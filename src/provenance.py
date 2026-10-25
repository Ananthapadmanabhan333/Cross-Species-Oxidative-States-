"""Run manifests: what was requested, what was received, what was written, with checksums."""

from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path
from typing import Any

from .data_collection.http_cache import sha256_file

TRACKED_PACKAGES = ("networkx", "pandas", "numpy", "PyYAML", "requests")


def software_versions() -> dict[str, str]:
    versions = {"python": sys.version.split()[0], "platform": platform.platform()}
    for package in TRACKED_PACKAGES:
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = "not installed"
    try:
        versions["extremonet"] = metadata.version("extremonet")
    except metadata.PackageNotFoundError:
        from . import __version__
        versions["extremonet"] = f"{__version__} (not installed)"
    return versions


def relative_posix(path: Path, root: Path) -> str:
    """Forward-slash path relative to `root` (portable across OSes), or absolute if outside it."""
    path, root = path.resolve(), root.resolve()
    return path.relative_to(root).as_posix() if path.is_relative_to(root) else path.as_posix()


def file_record(path: Path, root: Path) -> dict[str, Any]:
    return {
        "path": relative_posix(path, root),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def collect_outputs(directories: list[Path], root: Path, exclude: tuple[Path, ...] = ()) -> list[dict[str, Any]]:
    """Checksum every file under `directories`, except paths under `exclude` (e.g. a log that is
    still being written)."""
    records = []
    excluded = tuple(p.resolve() for p in exclude)
    for directory in directories:
        if directory.is_dir():
            for path in sorted(p for p in directory.rglob("*") if p.is_file()):
                if any(path.resolve().is_relative_to(e) for e in excluded):
                    continue
                records.append(file_record(path, root))
    return records


def write_manifest(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"manifest_written_at_utc": datetime.now(timezone.utc).isoformat(), **payload}
    path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")


def write_checksums(path: Path, records: list[dict[str, Any]]) -> None:
    """sha256sum-compatible listing (verify with `sha256sum -c` from the project root)."""
    lines = [f"{r['sha256']}  {r['path']}" for r in records]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
