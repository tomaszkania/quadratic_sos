"""Representability and exact row-wise square enumeration.

This module implements:
- Peters' representability criterion (Algorithm 1)
- exact row-wise enumeration of the search set ``B(alpha)``
"""

from __future__ import annotations

from typing import Dict, List, Optional, Set, Tuple

from quadratic_sos.ring import RealQuadraticOrder

Pair = Tuple[int, int]


def is_representable(order: RealQuadraticOrder, alpha: Pair) -> bool:
    """Decide whether ``alpha`` is a sum of squares in ``O_D``.

    This is Proposition 3.1 in the paper: Peters' criterion rewritten as a
    nearest-integer inequality, using exact integer arithmetic only.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    alpha : tuple[int, int]
        Coefficient pair representing ``alpha``.

    Returns
    -------
    bool
        ``True`` if and only if ``alpha`` is representable by integral squares.
    """
    if not order.is_totally_positive(alpha):
        return False

    a, b = alpha
    nrm = order.norm(alpha)

    if order.is_mod1:
        C = 2 * a + b
        c0 = _nearest_with_parity(C, order.D, b % 2)
        return (order.D * c0 - C) ** 2 <= 4 * nrm

    if b % 2 != 0:
        return False
    c0 = _nearest_integer_to_ratio(a, 2 * order.D)
    return (2 * order.D * c0 - a) ** 2 <= nrm


def _nearest_integer_to_ratio(num: int, den: int) -> int:
    """Return the integer nearest to ``num / den``.

    Parameters
    ----------
    num, den : int
        Integer numerator and positive denominator.

    Returns
    -------
    int
        An integer minimising ``|den*c - num|``.
    """
    q, r = divmod(num, den)
    if 2 * r >= den:
        return q + 1
    return q


def _nearest_with_parity(C: int, D: int, parity: int) -> int:
    """Return ``c ≡ parity (mod 2)`` minimising ``|D*c - C|``."""
    c_near = _nearest_integer_to_ratio(C, D)
    if c_near % 2 == parity:
        return c_near
    c_lo, c_hi = c_near - 1, c_near + 1
    if abs(D * c_lo - C) <= abs(D * c_hi - C):
        return c_lo
    return c_hi


def enumerate_search_set(order: RealQuadraticOrder, alpha: Pair) -> List[Pair]:
    """Enumerate ``B(alpha)`` exactly, row by row.

    The implementation avoids floating-point arithmetic completely. For each
    fixed row ``v``, it first computes an exact trace-based enclosure for the
    admissible integers ``u`` and then refines this to the exact row interval by
    directed binary search using only exact embedding comparisons.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    alpha : tuple[int, int]
        Totally positive target element.

    Returns
    -------
    list[tuple[int, int]]
        All coefficient pairs ``(u, v)`` such that
        ``(u + v*omega_D)^2 <= alpha`` in both real embeddings.
    """
    roots: List[Pair] = []
    v_max = order.trace_v_bound(alpha)
    for v in range(-v_max, v_max + 1):
        interval = _exact_row_interval(order, alpha, v)
        if interval is None:
            continue
        u_min, u_max = interval
        roots.extend((u, v) for u in range(u_min, u_max + 1))
    return roots


def relevant_squares(order: RealQuadraticOrder, alpha: Pair) -> List[Pair]:
    """Return the distinct squares arising from ``B(alpha)``.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    alpha : tuple[int, int]
        Target element.

    Returns
    -------
    list[tuple[int, int]]
        Sorted list ``S(alpha) = {x^2 : x in B(alpha)}``.
    """
    seen: Set[Pair] = set()
    result: List[Pair] = []
    for root in enumerate_search_set(order, alpha):
        sq = order.sqr(root)
        if sq not in seen:
            seen.add(sq)
            result.append(sq)
    result.sort()
    return result


def relevant_squares_with_roots(order: RealQuadraticOrder, alpha: Pair) -> Dict[Pair, Pair]:
    """Return a map ``x^2 -> x`` for the search set ``B(alpha)``.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    alpha : tuple[int, int]
        Target element.

    Returns
    -------
    dict[tuple[int, int], tuple[int, int]]
        Mapping from each distinct square in ``S(alpha)`` to one chosen root.
    """
    mapping: Dict[Pair, Pair] = {}
    for root in enumerate_search_set(order, alpha):
        sq = order.sqr(root)
        mapping.setdefault(sq, root)
    return mapping


def _square_le_both(order: RealQuadraticOrder, alpha: Pair, x: Pair) -> bool:
    """Return ``True`` when ``x^2 <= alpha`` in both embeddings."""
    sq = order.sqr(x)
    return order.sigma1_le(sq, alpha) and order.sigma2_le(sq, alpha)


def _row_direction(order: RealQuadraticOrder, alpha: Pair, x: Pair) -> int:
    """Describe where the admissible row interval lies relative to ``x``.

    Returns
    -------
    int
        ``0`` if ``x`` is admissible,
        ``1`` if the admissible interval lies strictly to the right,
        ``-1`` if it lies strictly to the left,
        ``2`` if the two embedding intervals are disjoint on this row.
    """
    sq = order.sqr(x)
    inside1 = order.sigma1_le(sq, alpha)
    inside2 = order.sigma2_le(sq, alpha)
    if inside1 and inside2:
        return 0

    dir1 = 0
    if not inside1:
        dir1 = 1 if order._sigma1_sign(x) < 0 else -1

    dir2 = 0
    if not inside2:
        dir2 = 1 if order._sigma2_sign(x) < 0 else -1

    if dir1 and dir2 and dir1 != dir2:
        return 2
    return dir1 or dir2


def _exact_row_interval(
    order: RealQuadraticOrder,
    alpha: Pair,
    v: int,
) -> Optional[Tuple[int, int]]:
    """Compute the exact admissible interval of ``u`` on the row ``v``.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    alpha : tuple[int, int]
        Totally positive target element.
    v : int
        Row index.

    Returns
    -------
    tuple[int, int] or None
        Exact interval ``[u_min, u_max]`` of admissible integers ``u`` on row
        ``v``. Returns ``None`` when the row contains no admissible points.
    """
    bounds = order.trace_row_bounds(alpha, v)
    if bounds is None:
        return None
    left, right = bounds
    if left > right:
        return None

    found = _find_admissible_u(order, alpha, v, left, right)
    if found is None:
        return None

    u_min = _left_boundary(order, alpha, v, left, found)
    u_max = _right_boundary(order, alpha, v, found, right)
    return (u_min, u_max)


def _find_admissible_u(
    order: RealQuadraticOrder,
    alpha: Pair,
    v: int,
    left: int,
    right: int,
) -> Optional[int]:
    """Find one admissible integer ``u`` on a fixed row by directed bisection."""
    lo, hi = left, right
    while lo <= hi:
        mid = (lo + hi) // 2
        direction = _row_direction(order, alpha, (mid, v))
        if direction == 0:
            return mid
        if direction == 1:
            lo = mid + 1
            continue
        if direction == -1:
            hi = mid - 1
            continue
        return None
    return None


def _left_boundary(
    order: RealQuadraticOrder,
    alpha: Pair,
    v: int,
    left: int,
    found: int,
) -> int:
    """Return the smallest admissible ``u`` on a fixed row."""
    lo, hi = left, found
    while lo < hi:
        mid = (lo + hi) // 2
        if _square_le_both(order, alpha, (mid, v)):
            hi = mid
        else:
            lo = mid + 1
    return lo


def _right_boundary(
    order: RealQuadraticOrder,
    alpha: Pair,
    v: int,
    found: int,
    right: int,
) -> int:
    """Return the largest admissible ``u`` on a fixed row."""
    lo, hi = found, right
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if _square_le_both(order, alpha, (mid, v)):
            lo = mid
        else:
            hi = mid - 1
    return lo
