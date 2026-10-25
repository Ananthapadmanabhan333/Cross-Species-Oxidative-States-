"""HTTP access with an on-disk, provenance-preserving cache.

Every response that reaches the pipeline is written to ``raw_dir`` exactly as received, next to a
``.meta.json`` sidecar holding the request (method, URL, parameters), retrieval time, HTTP status,
selected response headers and the SHA-256 of the body. Files are never overwritten: a refresh adds
a new timestamped copy, so earlier raw inputs stay available for audit.
"""

from __future__ import annotations

import hashlib
import json
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

import requests

from ..logging_utils import get_logger

log = get_logger("http")

RETRYABLE_STATUS = {429, 500, 502, 503, 504}
KEPT_HEADERS = (
    "content-type", "date", "x-total-results", "x-uniprot-release", "x-uniprot-release-date",
    "x-api-deployment-date",
)


class DataSourceError(RuntimeError):
    """The remote source answered, but with a non-retryable error (e.g. HTTP 400/404)."""


class DataSourceUnavailable(RuntimeError):
    """The remote source could not be reached after all retries."""


class OfflineCacheMiss(RuntimeError):
    """Offline mode was requested but the needed response has never been cached."""


class CacheIntegrityError(RuntimeError):
    """A cached file no longer matches the checksum recorded when it was downloaded."""


@dataclass
class RawResponse:
    text: str
    meta: dict[str, Any]
    data_path: Path
    meta_path: Path
    served_from: str  # "network" | "cache" | "stale_cache"
    headers: dict[str, str] = field(default_factory=dict)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def request_key(method: str, url: str, params: dict[str, Any]) -> str:
    """Deterministic identifier for a request; identical requests share cached responses."""
    payload = json.dumps(
        {"method": method.upper(), "url": url, "params": {k: str(v) for k, v in sorted(params.items())}},
        sort_keys=True,
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:20]


class CachedHttpClient:
    def __init__(
        self,
        raw_dir: Path,
        mode: str = "online",
        timeout_seconds: float = 60,
        max_retries: int = 4,
        backoff_seconds: float = 2.0,
        min_interval_seconds: float = 1.0,
        allow_stale_cache_on_failure: bool = False,
        user_agent: str = "extremonet",
        session: Any = None,
        sleep: Callable[[float], None] = time.sleep,
    ) -> None:
        if mode not in {"online", "refresh", "offline"}:
            raise ValueError(f"unknown mode {mode!r}")
        self.raw_dir = Path(raw_dir)
        self.mode = mode
        self.timeout = timeout_seconds
        self.max_retries = max_retries
        self.backoff = backoff_seconds
        self.min_interval = min_interval_seconds
        self.allow_stale = allow_stale_cache_on_failure
        self.session = session if session is not None else requests.Session()
        self.headers = {"User-Agent": user_agent}
        self._sleep = sleep
        self._last_request_at: float | None = None
        self.requests_log: list[dict[str, Any]] = []  # every response served, for the manifest

    @classmethod
    def from_config(cls, config: dict[str, Any], **kwargs: Any) -> "CachedHttpClient":
        http = config["http"]
        return cls(
            raw_dir=Path(config["paths"]["raw_dir"]),
            mode=config["mode"],
            timeout_seconds=http["timeout_seconds"],
            max_retries=http["max_retries"],
            backoff_seconds=http.get("backoff_seconds", 2.0),
            min_interval_seconds=http.get("min_interval_seconds", 1.0),
            allow_stale_cache_on_failure=http.get("allow_stale_cache_on_failure", False),
            user_agent=http.get("user_agent", "extremonet"),
            **kwargs,
        )

    # ------------------------------------------------------------------ public API
    def fetch(
        self,
        source: str,
        endpoint: str,
        url: str,
        params: dict[str, Any],
        method: str = "GET",
        ext: str = "txt",
        description: str = "",
        headers: dict[str, str] | None = None,
    ) -> RawResponse:
        key = request_key(method, url, params)
        cache_dir = self.raw_dir / source / endpoint
        cached = self._latest_cached(cache_dir, key)

        if self.mode == "offline":
            if cached is None:
                raise OfflineCacheMiss(
                    f"offline mode: no cached response for {source}/{endpoint} (key {key}). "
                    "Run once with mode=online to acquire it, or copy a previously acquired data/raw "
                    "directory into place."
                )
            return self._record(self._load(cached, "cache"))

        if self.mode == "online" and cached is not None:
            return self._record(self._load(cached, "cache"))

        try:
            response = self._request_with_retries(method, url, params, headers)
        except DataSourceUnavailable:
            if cached is not None and self.allow_stale:
                log.warning(
                    "source unavailable; using older cached response (allow_stale_cache_on_failure=true)",
                    extra={"ctx": {"source": source, "endpoint": endpoint, "key": key}},
                )
                return self._record(self._load(cached, "stale_cache"))
            raise

        return self._record(self._store(cache_dir, key, source, endpoint, method, url, params, ext,
                                        description, response))

    # ------------------------------------------------------------------ internals
    def _throttle(self) -> None:
        if self._last_request_at is not None and self.min_interval > 0:
            wait = self.min_interval - (time.monotonic() - self._last_request_at)
            if wait > 0:
                self._sleep(wait)
        self._last_request_at = time.monotonic()

    def _request_with_retries(self, method: str, url: str, params: dict[str, Any],
                              headers: dict[str, str] | None = None) -> Any:
        last_error = "unknown error"
        for attempt in range(self.max_retries + 1):
            self._throttle()
            try:
                kwargs: dict[str, Any] = {"timeout": self.timeout, "headers": {**self.headers, **(headers or {})}}
                if method.upper() == "POST":
                    kwargs["data"] = params
                else:
                    kwargs["params"] = params
                response = self.session.request(method.upper(), url, **kwargs)
            except (requests.ConnectionError, requests.Timeout) as exc:
                last_error = f"{type(exc).__name__}: {exc}"
            else:
                if response.status_code == 200:
                    return response
                body = (response.text or "")[:300]
                if response.status_code not in RETRYABLE_STATUS:
                    raise DataSourceError(f"{method} {url} returned HTTP {response.status_code}: {body}")
                last_error = f"HTTP {response.status_code}: {body}"
                retry_after = response.headers.get("Retry-After") if response.headers else None
                if retry_after and retry_after.isdigit() and attempt < self.max_retries:
                    log.warning("retrying after server-requested delay",
                                extra={"ctx": {"url": url, "attempt": attempt + 1, "wait_s": int(retry_after)}})
                    self._sleep(float(retry_after))
                    continue

            if attempt < self.max_retries:
                wait = self.backoff * (2 ** attempt)
                log.warning("request failed; retrying",
                            extra={"ctx": {"url": url, "attempt": attempt + 1, "error": last_error, "wait_s": wait}})
                self._sleep(wait)

        raise DataSourceUnavailable(
            f"{method} {url} failed after {self.max_retries + 1} attempts ({last_error}). Recovery: "
            "re-run later, or run with mode=offline using a previously acquired data/raw directory."
        )

    def _store(self, cache_dir: Path, key: str, source: str, endpoint: str, method: str, url: str,
               params: dict[str, Any], ext: str, description: str, response: Any) -> RawResponse:
        cache_dir.mkdir(parents=True, exist_ok=True)
        retrieved_at = datetime.now(timezone.utc)
        stamp = retrieved_at.strftime("%Y%m%dT%H%M%S%fZ")
        body = response.content if isinstance(response.content, bytes) else str(response.text).encode("utf-8")
        # The OS clock can repeat a timestamp for fetches milliseconds apart; never overwrite an
        # earlier raw file. "_NN" sorts after "." so the newest copy is still found last by name.
        n = 1
        while (cache_dir / f"{key}__{stamp}.meta.json").exists():
            n += 1
            stamp = f"{retrieved_at.strftime('%Y%m%dT%H%M%S%fZ')}_{n:02d}"
        data_path = cache_dir / f"{key}__{stamp}.{ext}"
        meta_path = cache_dir / f"{key}__{stamp}.meta.json"
        data_path.write_bytes(body)
        headers = {k: v for k, v in ((h, response.headers.get(h)) for h in KEPT_HEADERS)
                   if v is not None} if response.headers else {}
        meta = {
            "source": source,
            "endpoint": endpoint,
            "description": description,
            "request_key": key,
            "method": method.upper(),
            "url": url,
            "params": {k: str(v) for k, v in sorted(params.items())},
            "retrieved_at_utc": retrieved_at.isoformat(),
            "http_status": response.status_code,
            "response_headers": headers,
            "bytes": len(body),
            "sha256": sha256_bytes(body),
            "data_file": data_path.name,
        }
        meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
        log.info("downloaded", extra={"ctx": {"source": source, "endpoint": endpoint, "bytes": len(body),
                                              "file": str(data_path)}})
        return RawResponse(body.decode("utf-8"), meta, data_path, meta_path, "network", headers)

    @staticmethod
    def _latest_cached(cache_dir: Path, key: str) -> Path | None:
        if not cache_dir.is_dir():
            return None
        metas = sorted(cache_dir.glob(f"{key}__*.meta.json"))
        return metas[-1] if metas else None

    @staticmethod
    def _load(meta_path: Path, served_from: str) -> RawResponse:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        data_path = meta_path.parent / meta["data_file"]
        if not data_path.is_file():
            raise CacheIntegrityError(f"cached data file missing: {data_path}")
        actual = sha256_file(data_path)
        if actual != meta["sha256"]:
            raise CacheIntegrityError(
                f"checksum mismatch for {data_path}: recorded {meta['sha256']}, found {actual}"
            )
        log.info("loaded from cache", extra={"ctx": {"source": meta["source"], "endpoint": meta["endpoint"],
                                                     "retrieved_at": meta["retrieved_at_utc"]}})
        return RawResponse(data_path.read_bytes().decode("utf-8"), meta, data_path, meta_path, served_from,
                           meta.get("response_headers", {}))

    def _record(self, raw: RawResponse) -> RawResponse:
        self.requests_log.append({
            "source": raw.meta["source"],
            "endpoint": raw.meta["endpoint"],
            "description": raw.meta.get("description", ""),
            "served_from": raw.served_from,
            "url": raw.meta["url"],
            "method": raw.meta["method"],
            "params": raw.meta["params"],
            "retrieved_at_utc": raw.meta["retrieved_at_utc"],
            "sha256": raw.meta["sha256"],
            "data_file": str(raw.data_path),
        })
        return raw
