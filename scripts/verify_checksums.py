"""Verify a run's outputs and every cached raw API response against their recorded SHA-256.

Usage (from the project root):
    python scripts/verify_checksums.py                 # latest run
    python scripts/verify_checksums.py <run_id>
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(argv: list[str]) -> int:
    root = Path.cwd()
    run_id = argv[0] if argv else (root / "results" / "LATEST_RUN.txt").read_text().strip()
    failures = 0

    sums = root / "results" / run_id / "SHA256SUMS.txt"
    lines = [line for line in sums.read_text(encoding="utf-8").splitlines() if line]
    for line in lines:
        digest, rel = line.split("  ", 1)
        path = root / rel
        if not path.is_file() or sha256(path) != digest:
            print(f"FAILED output {rel}")
            failures += 1

    metas = sorted((root / "data" / "raw").rglob("*.meta.json"))
    for meta_path in metas:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        data = meta_path.parent / meta["data_file"]
        if not data.is_file() or sha256(data) != meta["sha256"]:
            print(f"FAILED raw {data}")
            failures += 1

    print(f"run {run_id}: {len(lines)} outputs, {len(metas)} raw responses checked, {failures} failures")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
