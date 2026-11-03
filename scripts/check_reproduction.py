"""Check that a re-run reproduces a reference run.

Every .tsv and .graphml table must be byte-identical. JSON summaries must agree to 12 significant
digits (floats are written at that precision; older runs stored full precision), ignoring the
fields that only identify the run. Files that record when or how a run was made (manifest,
checksums, logs, config copy, generated report and figures) are skipped.

Usage (from the project root):
    python scripts/check_reproduction.py <reference_run_dir> <new_run_dir>
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

SKIP_NAMES = {"manifest.json", "SHA256SUMS.txt", "config_used.yaml", "technical_report.md"}
SKIP_DIRS = {"logs", "figures"}
TABLE_SUFFIXES = {".tsv", ".graphml"}
RUN_IDENTITY_KEYS = {"run_id"}


def result_files(run: Path) -> dict[str, Path]:
    found = {}
    for path in sorted(run.rglob("*")):
        rel = path.relative_to(run)
        if (path.is_file() and path.name not in SKIP_NAMES and not SKIP_DIRS & set(rel.parts)
                and path.suffix in TABLE_SUFFIXES | {".json"}):
            found[rel.as_posix()] = path
    return found


def same_json(a: Any, b: Any) -> bool:
    if isinstance(a, dict) and isinstance(b, dict):
        keys = (set(a) | set(b)) - RUN_IDENTITY_KEYS
        return all(k in a and k in b and same_json(a[k], b[k]) for k in keys)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(same_json(x, y) for x, y in zip(a, b))
    if isinstance(a, float) or isinstance(b, float):
        if not isinstance(a, (int, float)) or not isinstance(b, (int, float)):
            return False
        if math.isnan(a) or math.isnan(b):
            return math.isnan(a) and math.isnan(b)
        return float(f"{a:.12g}") == float(f"{b:.12g}")
    return a == b


def same(path_a: Path, path_b: Path) -> bool:
    if path_a.suffix == ".json":
        return same_json(json.loads(path_a.read_text(encoding="utf-8")),
                         json.loads(path_b.read_text(encoding="utf-8")))
    return hashlib.sha256(path_a.read_bytes()).digest() == hashlib.sha256(path_b.read_bytes()).digest()


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__)
        return 2
    ref, new = (Path(a) for a in argv)
    ref_files, new_files = result_files(ref), result_files(new)
    missing = sorted(set(ref_files) ^ set(new_files))
    for rel in missing:
        print(f"MISSING in {'new' if rel in ref_files else 'reference'} run: {rel}")
    differ = [rel for rel in sorted(set(ref_files) & set(new_files)) if not same(ref_files[rel], new_files[rel])]
    for rel in differ:
        print(f"DIFFERS: {rel}")
    print(f"{ref.name} vs {new.name}: {len(ref_files)} result files, "
          f"{len(missing)} missing, {len(differ)} differing")
    return 1 if missing or differ else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
