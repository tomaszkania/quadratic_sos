#!/usr/bin/env python3
"""Brute-force validation of exact lengths on bounded trace ranges."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from quadratic_sos.experiments import validation_sweep


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-D", type=int, default=41, help="validate all squarefree D up to this bound")
    parser.add_argument("--trace-bound", type=int, default=24, help="maximum trace for exhaustive brute-force comparison")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary = validation_sweep(args.max_D, args.trace_bound)
    print(f"Squarefree fields checked : {summary.field_count}")
    print(f"Trace bound             : {summary.trace_bound}")
    print(f"Elements compared       : {summary.element_count}")
    print(f"Mismatches              : {summary.mismatches}")


if __name__ == "__main__":
    main()
