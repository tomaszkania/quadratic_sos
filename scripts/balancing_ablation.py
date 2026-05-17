#!/usr/bin/env python3
"""Balancing ablation on a Pell orbit.

This script compares the raw and balanced search geometry for elements
``alpha_k = u^(2k) * beta`` in a fixed real quadratic ring, where ``u`` is the
positive norm-one Pell unit used by the package.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from quadratic_sos import RealQuadraticOrder
from quadratic_sos.experiments import pell_orbit, search_profile
from quadratic_sos.enumeration import enumerate_search_set


def _median_time_ms(order: RealQuadraticOrder, alpha: tuple[int, int], repeats: int = 3) -> float:
    samples = []
    for _ in range(repeats):
        start = time.perf_counter()
        enumerate_search_set(order, alpha)
        samples.append((time.perf_counter() - start) * 1000.0)
    samples.sort()
    return samples[len(samples) // 2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--D", type=int, default=10, help="squarefree radicand")
    parser.add_argument("--a", type=int, default=18, help="rational coefficient of beta")
    parser.add_argument("--b", type=int, default=2, help="omega-coefficient of beta")
    parser.add_argument("--max-k", type=int, default=3, help="largest Pell-orbit index")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    order = RealQuadraticOrder(args.D)
    beta = (args.a, args.b)
    print(f"Balancing ablation in O_{order.D} on the orbit of {order.elem(*beta)!r}")
    print()
    print("| k | Tr(raw) | V_tr(raw) | rows(raw) | |B|(raw) | raw ms | Tr(bal) | V_tr(bal) | rows(bal) | |B|(bal) | bal ms |")
    print("| -: | ------: | --------: | --------: | -------: | -----: | ------: | --------: | --------: | -------: | -----: |")
    for k, alpha in enumerate(pell_orbit(order, beta, args.max_k)):
        raw = search_profile(order, alpha, balance=False)
        bal = search_profile(order, alpha, balance=True)
        raw_ms = _median_time_ms(order, raw.search_alpha)
        bal_ms = _median_time_ms(order, bal.search_alpha)
        print(
            f"| {k} | {raw.trace} | {raw.trace_v_bound} | {raw.scanned_rows} | {raw.box_size} | {raw_ms:.3f} | "
            f"{bal.trace} | {bal.trace_v_bound} | {bal.scanned_rows} | {bal.box_size} | {bal_ms:.3f} |"
        )


if __name__ == "__main__":
    main()
