"""Utilities for computational experiments appearing in the paper.

The core package focuses on exact length and decomposition.  This module adds
small helpers used by the computational section: search-profile statistics,
Pell-orbit generation, and brute-force validation on small inputs.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
import math
from typing import Dict, List, Optional, Tuple

from quadratic_sos.enumeration import _exact_row_interval, enumerate_search_set, relevant_squares
from quadratic_sos.length import exact_length
from quadratic_sos.ring import RealQuadraticOrder
from quadratic_sos.tables import enumerate_totally_positive

Pair = Tuple[int, int]


@dataclass(frozen=True)
class SearchProfile:
    """Summary statistics for one exact search.

    Attributes
    ----------
    input_alpha : tuple[int, int]
        Original input element.
    search_alpha : tuple[int, int]
        Element actually searched; this equals ``input_alpha`` when balancing is
        disabled and a balanced unit-square multiple otherwise.
    balanced : bool
        Whether the balancing step was applied.
    trace : int
        Trace of ``search_alpha``.
    trace_v_bound : int
        Exact trace-based row bound ``V_tr``.
    scanned_rows : int
        Number of rows scanned, i.e. ``2*trace_v_bound + 1``.
    nonempty_rows : int
        Number of rows whose exact admissible interval is nonempty.
    box_size : int
        Cardinality ``|B(search_alpha)|`` of admissible roots.
    square_count : int
        Cardinality ``|S(search_alpha)|`` of distinct relevant squares.
    pair_sum_count : int
        Cardinality ``|T_2(search_alpha)|`` of distinct unordered pair sums.
    sigma_row : float
        The quantity ``Σ_row = sum_v log_2(1 + |I_v^tr|)`` over all trace rows.
    """

    input_alpha: Pair
    search_alpha: Pair
    balanced: bool
    trace: int
    trace_v_bound: int
    scanned_rows: int
    nonempty_rows: int
    box_size: int
    square_count: int
    pair_sum_count: int
    sigma_row: float

    def as_dict(self) -> Dict[str, object]:
        """Return the profile as a plain dictionary."""
        return asdict(self)


def pell_orbit(order: RealQuadraticOrder, beta: Pair, max_k: int) -> List[Pair]:
    """Return ``[beta, u^2 beta, u^4 beta, ..., u^{2 max_k} beta]``.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    beta : tuple[int, int]
        Base element.
    max_k : int
        Largest exponent index.

    Returns
    -------
    list[tuple[int, int]]
        Pell orbit of ``beta`` under even powers of the positive norm-one Pell
        unit returned by :attr:`RealQuadraticOrder.pell_unit`.
    """
    if max_k < 0:
        raise ValueError("max_k must be non-negative")
    unit_sq = order.sqr(order.pell_unit)
    out = [beta]
    current = beta
    for _ in range(max_k):
        current = order.mul(current, unit_sq)
        out.append(current)
    return out


def search_profile(order: RealQuadraticOrder, alpha: Pair, *, balance: bool = True) -> SearchProfile:
    """Compute the search statistics used in the computational section.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    alpha : tuple[int, int]
        Input element.
    balance : bool, default=True
        Whether to replace ``alpha`` by a balanced unit-square multiple before
        collecting statistics.

    Returns
    -------
    SearchProfile
        Exact structural statistics for the search performed on ``alpha``.
    """
    search_alpha, _ = order.balance_by_pell_unit(alpha) if balance else (alpha, order.one)
    v_max = order.trace_v_bound(search_alpha)
    scanned_rows = 2 * v_max + 1
    nonempty_rows = 0
    box_size = 0
    sigma_row = 0.0

    for v in range(-v_max, v_max + 1):
        bounds = order.trace_row_bounds(search_alpha, v)
        if bounds is None:
            continue
        left, right = bounds
        width = max(0, right - left + 1)
        if width:
            sigma_row += math.log2(1 + width)
        row = _exact_row_interval(order, search_alpha, v)
        if row is None:
            continue
        nonempty_rows += 1
        box_size += row[1] - row[0] + 1

    squares = relevant_squares(order, search_alpha)
    # Distinct pair sums are more informative than raw unordered pairs.
    t2 = set()
    for idx, s1 in enumerate(squares):
        for s2 in squares[idx:]:
            t2.add(order.add(s1, s2))

    return SearchProfile(
        input_alpha=alpha,
        search_alpha=search_alpha,
        balanced=balance,
        trace=order.trace(search_alpha),
        trace_v_bound=v_max,
        scanned_rows=scanned_rows,
        nonempty_rows=nonempty_rows,
        box_size=box_size,
        square_count=len(squares),
        pair_sum_count=len(t2),
        sigma_row=sigma_row,
    )


def brute_force_relevant_squares(order: RealQuadraticOrder, alpha: Pair) -> List[Pair]:
    """Enumerate ``S(alpha)`` by a deliberately naïve coefficient-box scan.

    This routine is intentionally independent of the exact row-search engine and
    is used only for validation on small trace ranges.
    """
    tr = order.trace(alpha)
    u_max = tr + 1
    v_max = tr + 1
    out = set()
    for u in range(-u_max, u_max + 1):
        for v in range(-v_max, v_max + 1):
            x = (u, v)
            sq = order.sqr(x)
            if order.sigma1_le(sq, alpha) and order.sigma2_le(sq, alpha):
                out.add(sq)
    return sorted(out)


def brute_force_exact_length(order: RealQuadraticOrder, alpha: Pair) -> Optional[int]:
    """Compute the exact length by an intentionally naive validator.

    The routine is designed only for validation on modest trace ranges. It scans
    a coefficient box and tests all sums of at most five relevant squares. It
    deliberately avoids Peters' representability criterion, the exact row-search
    engine, and the closed-form Pythagoras-number classification.
    """
    squares = brute_force_relevant_squares(order, alpha)
    square_set = set(squares)

    if alpha in square_set:
        return 1

    for square in squares:
        if order.sub(alpha, square) in square_set:
            return 2

    pair_sums = set()
    for index, first in enumerate(squares):
        for second in squares[index:]:
            pair_sums.add(order.add(first, second))

    for square in squares:
        if order.sub(alpha, square) in pair_sums:
            return 3

    for pair_sum in pair_sums:
        if order.sub(alpha, pair_sum) in pair_sums:
            return 4

    for square in squares:
        remainder = order.sub(alpha, square)
        for pair_sum in pair_sums:
            if order.sub(remainder, pair_sum) in pair_sums:
                return 5

    return None


@dataclass(frozen=True)
class ValidationSummary:
    """Summary of a brute-force validation sweep."""

    field_count: int
    trace_bound: int
    element_count: int
    mismatches: int


def validation_sweep(max_D: int, trace_bound: int) -> ValidationSummary:
    """Validate the exact-length routine against brute force on small inputs.

    Parameters
    ----------
    max_D : int
        Validate all squarefree ``D`` with ``2 <= D <= max_D``.
    trace_bound : int
        Compare all nonzero totally positive elements of trace at most this
        bound.

    Returns
    -------
    ValidationSummary
        Counts of fields, elements, and mismatches.
    """
    Ds = [D for D in range(2, max_D + 1) if _is_squarefree_non_square(D)]
    elements = 0
    mismatches = 0
    for D in Ds:
        order = RealQuadraticOrder(D)
        for alpha in enumerate_totally_positive(order, trace_bound):
            elements += 1
            if exact_length(order, alpha) != brute_force_exact_length(order, alpha):
                mismatches += 1
    return ValidationSummary(
        field_count=len(Ds),
        trace_bound=trace_bound,
        element_count=elements,
        mismatches=mismatches,
    )


def _is_squarefree_non_square(D: int) -> bool:
    """Return ``True`` when ``D`` is squarefree and not a square."""
    s = math.isqrt(D)
    if s * s == D:
        return False
    n = D
    p = 2
    while p * p <= n:
        exp = 0
        while n % p == 0:
            n //= p
            exp += 1
        if exp >= 2:
            return False
        p += 1
    return True
