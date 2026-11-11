"""Restore the KEGG responses that the repository does not redistribute.

KEGG's terms do not allow its data to be redistributed, so the repository keeps only the
provenance record (`*.meta.json`: URL, retrieval time, size and SHA-256) of each KEGG response,
not the response itself. This script downloads each response again from the KEGG REST API
(https://rest.kegg.jp) and stores it under the file name the record expects, but only if its
SHA-256 matches the record, i.e. if KEGG still serves exactly the data the reference runs used.

Exit status: 0 all responses present and verified; 3 KEGG now serves different data (the
reference runs cannot be reproduced byte for byte; use --accept-current to cache the current
data as a new response instead); 1 any other error.

    python scripts/fetch_kegg.py [--accept-current]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
KEGG_DIR = ROOT / "data" / "raw" / "kegg"
USER_AGENT = "extremonet-kegg-fetch (research pipeline)"


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def download(url: str) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=120) as response:  # noqa: S310 (fixed https URL)
        if response.status != 200:
            raise RuntimeError(f"{url} returned HTTP {response.status}")
        return response.read()


def store_current(meta_path: Path, meta: dict, body: bytes) -> Path:
    """Cache the current KEGG data as a new response, named like the pipeline's own cache files."""
    now = datetime.now(timezone.utc)
    stamp = now.strftime("%Y%m%dT%H%M%S%fZ")
    key = meta["request_key"]
    data_path = meta_path.parent / f"{key}__{stamp}.tsv"
    data_path.write_bytes(body)
    new_meta = dict(meta, retrieved_at_utc=now.isoformat(), bytes=len(body), sha256=sha256(body),
                    data_file=data_path.name, response_headers={})
    (meta_path.parent / f"{key}__{stamp}.meta.json").write_text(json.dumps(new_meta, indent=2),
                                                                 encoding="utf-8")
    return data_path


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--accept-current", action="store_true",
                        help="if KEGG's data changed, cache the current data as a new response")
    args = parser.parse_args(argv)

    metas = sorted(KEGG_DIR.rglob("*.meta.json"))
    if not metas:
        print(f"no KEGG provenance records under {KEGG_DIR}")
        return 1
    changed = 0
    for meta_path in metas:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        data_path = meta_path.parent / meta["data_file"]
        if data_path.is_file() and sha256(data_path.read_bytes()) == meta["sha256"]:
            print(f"present   {data_path.relative_to(ROOT)}")
            continue
        body = download(meta["url"])
        time.sleep(1)  # KEGG asks clients not to send rapid bursts of requests
        if sha256(body) == meta["sha256"]:
            data_path.write_bytes(body)
            print(f"restored  {data_path.relative_to(ROOT)} (matches the recorded SHA-256)")
        elif args.accept_current:
            new_path = store_current(meta_path, meta, body)
            print(f"changed   {meta['url']}: cached the current data as {new_path.relative_to(ROOT)}")
        else:
            changed += 1
            print(f"CHANGED   {meta['url']}: KEGG now serves different data "
                  f"(recorded {meta['sha256'][:12]}, now {sha256(body)[:12]}); nothing stored")
    if changed:
        print("KEGG data differ from the reference snapshot, so the reference runs cannot be reproduced "
              "byte for byte. Re-run with --accept-current to analyse the current KEGG release.")
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
