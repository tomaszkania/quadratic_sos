#!/usr/bin/env python3
"""Cross-field exact-length distributions for bounded trace."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from quadratic_sos import RealQuadraticOrder
from quadratic_sos.tables import length_distribution


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fields", type=int, nargs="+", default=[5, 10, 13, 29, 53, 101], help="list of squarefree fields")
    parser.add_argument("--trace-bound", type=int, default=80, help="maximum trace")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    print(f"Cross-field distributions for trace <= {args.trace_bound}")
    print()
    header = "| D | total | #1 | #2 | #3 | #4 | #5 | #∞ |"
    rule = "| -: | ----: | -: | -: | -: | -: | -: | -: |"
    print(header)
    print(rule)
    for D in args.fields:
        order = RealQuadraticOrder(D)
        counts, _ = length_distribution(order, args.trace_bound)
        total = sum(counts.values())
        print(
            f"| {D} | {total} | {counts.get(1, 0)} | {counts.get(2, 0)} | {counts.get(3, 0)} | "
            f"{counts.get(4, 0)} | {counts.get(5, 0)} | {counts.get(None, 0)} |"
        )


if __name__ == "__main__":
    main()
