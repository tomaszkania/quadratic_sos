"""Cross-validate batched and naive bounded truant routines on a small family."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from quadratic_diagonal import RealQuadraticOrder, bounded_truants_batched, bounded_truants_naive


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fields', nargs='*', type=int, default=[5, 6, 10, 13, 21])
    parser.add_argument('--trace-bounds', nargs='*', type=int, default=[8, 12, 16])
    args = parser.parse_args()

    total = 0
    for D in args.fields:
        order = RealQuadraticOrder(D)
        for coeffs in ([(1, 0), (1, 0)], [(1, 0), (1, 0), (1, 0)]):
            for B in args.trace_bounds:
                batched = bounded_truants_batched(order, coeffs, B)
                naive = bounded_truants_naive(order, coeffs, B, method='mitm')
                if batched != naive:
                    raise AssertionError(f'Mismatch for D={D}, coeffs={coeffs}, B={B}')
                total += 1
    print(f'Validated {total} bounded-truant instances successfully.')


if __name__ == '__main__':
    main()
