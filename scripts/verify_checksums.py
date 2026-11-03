"""Verify run outputs and every cached raw API response against their recorded SHA-256.

Usage (from the project root):
    python scripts/verify_checksums.py                 # latest Phase 1 run
    python scripts/verify_checksums.py <run_id>        # one Phase 1 run
    python scripts/verify_checksums.py --all           # every run under results/ (used by CI)
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sums(root: Path, sums: Path) -> tuple[int, int]:
    """Check one SHA256SUMS.txt; paths in it are relative to the project root."""
    lines = [line for line in sums.read_text(encoding="utf-8").splitlines() if line]
    failures = 0
    for line in lines:
        digest, rel = line.split("  ", 1)
        path = root / rel
        if not path.is_file() or sha256(path) != digest:
            print(f"FAILED output {rel}")
            failures += 1
    return len(lines), failures


def check_raw(root: Path) -> tuple[int, int]:
    metas = sorted((root / "data" / "raw").rglob("*.meta.json"))
    failures = 0
    for meta_path in metas:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        data = meta_path.parent / meta["data_file"]
        if not data.is_file() or sha256(data) != meta["sha256"]:
            print(f"FAILED raw {data}")
            failures += 1
    return len(metas), failures


def main(argv: list[str]) -> int:
    root = Path.cwd()
    if argv == ["--all"]:
        sums_files = sorted((root / "results").rglob("SHA256SUMS.txt"))
    else:
        run_id = argv[0] if argv else (root / "results" / "LATEST_RUN.txt").read_text().strip()
        sums_files = [root / "results" / run_id / "SHA256SUMS.txt"]

    failures = 0
    for sums in sums_files:
        n, failed = check_sums(root, sums)
        failures += failed
        print(f"run {sums.parent.relative_to(root).as_posix()}: {n} outputs, {failed} failures")
    n_raw, failed = check_raw(root)
    failures += failed
    print(f"raw cache: {n_raw} responses, {failed} failures")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
