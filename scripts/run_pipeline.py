"""Convenience wrapper; equivalent to `python -m extremonet`. Accepts the same arguments."""

from extremonet.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
