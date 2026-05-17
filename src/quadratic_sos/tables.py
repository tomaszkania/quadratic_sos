"""Bulk enumeration: compute exact-length distributions over O_D."""

from __future__ import annotations

from collections import Counter
from typing import Dict, List, Optional, Tuple

from quadratic_sos.ring import RealQuadraticOrder
from quadratic_sos.length import exact_length

Pair = Tuple[int, int]


def enumerate_totally_positive(order: RealQuadraticOrder, trace_bound: int) -> List[Pair]:
    """All totally positive alpha ∈ O_D with 1 ≤ Tr(alpha) ≤ trace_bound."""
    elems: List[Pair] = []
    D = order.D

    if order.is_mod1:
        # Tr(a + b*ω) = 2a + b.  Need 1 ≤ 2a+b ≤ trace_bound, tp.
        for tr in range(1, trace_bound + 1):
            # 2a + b = tr ⟹ b = tr - 2a
            # Need b to range so that alpha is tp.
            # a can be anything with b = tr - 2a integer.
            # sigma_1 > 0 and sigma_2 > 0 checked below.
            for a in range(-(trace_bound), trace_bound + 1):
                b = tr - 2 * a
                x = (a, b)
                if order.is_totally_positive(x):
                    elems.append(x)
    else:
        # Tr(a + b√D) = 2a.  Need 1 ≤ 2a ≤ trace_bound, so a ≥ 1.
        max_a = trace_bound // 2
        for a in range(1, max_a + 1):
            # |b| bounded by total positivity: a + b√D > 0 and a - b√D > 0
            # ⟹ |b| < a/√D.  Use a generous bound.
            b_bound = a  # always safe since √D > 1
            for b in range(-b_bound, b_bound + 1):
                x = (a, b)
                if order.is_totally_positive(x):
                    elems.append(x)

    elems.sort(key=lambda x: (order.trace(x), x[0], x[1]))
    return elems


def length_distribution(
    order: RealQuadraticOrder, trace_bound: int
) -> Tuple[Counter, Dict[Optional[int], Pair]]:
    """Compute length distribution and first examples.

    Returns
    -------
    counts : Counter
        Maps length (1–5 or None for ∞) to count.
    first : dict
        Maps length to the first element (by trace order) achieving it.
    """
    elems = enumerate_totally_positive(order, trace_bound)
    counts: Counter = Counter()
    first: Dict[Optional[int], Pair] = {}

    for alpha in elems:
        ell = exact_length(order, alpha)
        counts[ell] += 1
        if ell not in first:
            first[ell] = alpha

    return counts, first
