#!/usr/bin/env python3
"""Reproduce the paper's exact-length table using the package implementation."""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path
from typing import Optional

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from quadratic_sos import RealQuadraticOrder, exact_length
from quadratic_sos.tables import enumerate_totally_positive


Pair = tuple[int, int]


def format_pair(order: RealQuadraticOrder, pair: Pair) -> str:
    """Return a readable string for an element given by a coefficient pair."""
    return repr(order.elem(*pair))


def build_table(order: RealQuadraticOrder, trace_bound: int) -> None:
    """Print the exact-length distribution up to the given trace bound."""
    elems = enumerate_totally_positive(order, trace_bound)
    counts: Counter[Optional[int]] = Counter()
    first: dict[Optional[int], Pair] = {}

    for alpha in elems:
        ell = exact_length(order, alpha)
        counts[ell] += 1
        first.setdefault(ell, alpha)

    print(f"D = {order.D}, trace bound = {trace_bound}")
    print(f"Total totally positive elements counted: {len(elems)}")
    print()
    print(f"{'length':>8} | {'count':>5} | smallest trace-lex example")
    print("-" * 44)
    for key in [1, 2, 3, 4, 5, None]:
        label = "∞" if key is None else str(key)
        example = format_pair(order, first[key]) if key in first else '--'
        print(f"{label:>8} | {counts.get(key, 0):>5} | {example}")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--D", type=int, default=10, help="squarefree positive integer D")
    parser.add_argument("--trace-bound", type=int, default=40, help="upper bound for the trace")
    return parser.parse_args()


def main() -> None:
    """Run the table computation."""
    args = parse_args()
    build_table(RealQuadraticOrder(args.D), args.trace_bound)


if __name__ == "__main__":
    main()
