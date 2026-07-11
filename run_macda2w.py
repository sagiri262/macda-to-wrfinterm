#!/usr/bin/env python3
"""Command line entry point for MACDA v2.0 to WRF intermediate conversion."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from macda2wrf.config import load_config
from macda2wrf.converter import MacdaConverter


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-c",
        "--config",
        default="conf/config.MACDA-v2.ini",
        help="Path to the MACDA conversion config file.",
    )
    parser.add_argument(
        "--time-index",
        type=int,
        action="append",
        help="Process one time index. Repeat to process several explicit indices.",
    )
    parser.add_argument(
        "--max-times",
        type=int,
        help="Override config max_times for this run.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Open the file, print planned fields/times, and do not write output.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        cfg = load_config(Path(args.config))
        converter = MacdaConverter(cfg)
        written = converter.run(
            time_indices=args.time_index,
            max_times_override=args.max_times,
            dry_run=args.dry_run,
        )
    except Exception as exc:
        print(f"[macda2w] error: {exc}", file=sys.stderr)
        sys.exit(1)
    if written:
        print("[macda2w] wrote:")
        for path in written:
            print(f"  {path}")


if __name__ == "__main__":
    main()
