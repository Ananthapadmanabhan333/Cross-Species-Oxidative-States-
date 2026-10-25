"""Command-line entry point: `extremonet --config config/phase1_oxidative_stress.yaml [--mode offline]`."""

from __future__ import annotations

import argparse
import sys

from .config import ConfigError, load_config
from .data_collection.http_cache import (
    CacheIntegrityError, DataSourceError, DataSourceUnavailable, OfflineCacheMiss,
)
from .logging_utils import get_logger


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="extremonet", description=__doc__)
    parser.add_argument("--config", default="config/phase1_oxidative_stress.yaml", help="YAML config file")
    parser.add_argument("--mode", choices=["online", "refresh", "offline"], help="override config 'mode'")
    parser.add_argument("--threshold", type=float, help="override network.edge_score_threshold (0-1)")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    overrides = {}
    if args.mode:
        overrides["mode"] = args.mode
    if args.threshold is not None:
        overrides["network.edge_score_threshold"] = args.threshold
    try:
        config = load_config(args.config, overrides)
    except (ConfigError, FileNotFoundError) as exc:
        print(f"configuration error: {exc}", file=sys.stderr)
        return 2

    from .pipeline import run_pipeline  # imported late so config errors report fast

    try:
        result = run_pipeline(config)
    except (OfflineCacheMiss, DataSourceUnavailable, DataSourceError, CacheIntegrityError) as exc:
        get_logger("cli").error("data acquisition failed", extra={"ctx": {"error": str(exc)}})
        print(f"\nData acquisition failed: {exc}", file=sys.stderr)
        return 3
    print(f"\nRun {result['run_id']} complete.\n  results:  {result['results_dir']}\n"
          f"  manifest: {result['manifest']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
