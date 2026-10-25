"""Phase 2 entry point: `python -m extremonet.analysis --config config/phase2_analysis.yaml`."""

from __future__ import annotations

import argparse
import sys

from ..data_collection.http_cache import (
    CacheIntegrityError, DataSourceError, DataSourceUnavailable, OfflineCacheMiss,
)
from ..logging_utils import get_logger
from .inputs import InputIntegrityError


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="extremonet-analysis", description=__doc__)
    parser.add_argument("--config", default="config/phase2_analysis.yaml")
    parser.add_argument("--mode", choices=["online", "refresh", "offline"])
    parser.add_argument("--seed", type=int, help="override random_seed")
    args = parser.parse_args(argv)
    overrides = {}
    if args.mode:
        overrides["mode"] = args.mode
    if args.seed is not None:
        overrides["random_seed"] = args.seed

    from .phase2 import load_phase2_config, run_phase2

    try:
        cfg = load_phase2_config(args.config, overrides)
    except (ValueError, FileNotFoundError, KeyError) as exc:
        print(f"configuration error: {exc}", file=sys.stderr)
        return 2
    try:
        result = run_phase2(cfg)
    except InputIntegrityError as exc:
        print(f"Phase 1 inputs failed verification: {exc}", file=sys.stderr)
        return 4
    except (OfflineCacheMiss, DataSourceUnavailable, DataSourceError, CacheIntegrityError) as exc:
        get_logger("cli").error("data acquisition failed", extra={"ctx": {"error": str(exc)}})
        print(f"\nData acquisition failed: {exc}", file=sys.stderr)
        return 3
    print(f"\nPhase 2 run {result['run_id']} complete.\n  results:  {result['results_dir']}\n"
          f"  report:   {result['results_dir'] / 'technical_report.md'}\n  manifest: {result['manifest']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
