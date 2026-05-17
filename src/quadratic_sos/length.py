"""Exact integral length computation and decomposition recovery.

The main function :func:`exact_length` computes
``ell_D(alpha) in {1, 2, 3, 4, 5, None}``, where ``None`` represents
non-representability. The function :func:`decomposition` additionally returns
an explicit minimal sum-of-squares decomposition.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Set, Tuple
import warnings

from quadratic_sos.enumeration import is_representable, relevant_squares, relevant_squares_with_roots
from quadratic_sos.ring import RealQuadraticOrder

Pair = Tuple[int, int]
ZERO: Pair = (0, 0)


def exact_length(order: RealQuadraticOrder, alpha: Pair, *, balance: bool = True) -> Optional[int]:
    """Compute the exact integral length ``ell_D(alpha)`` for nonzero ``alpha``.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    alpha : tuple[int, int]
        Nonzero coefficient pair ``(a, b)`` representing ``a + b*omega_D``.

    Returns
    -------
    int or None
        The minimal number of squares, or ``None`` when ``alpha`` is not
        representable. When ``balance`` is ``True`` (the default), the search
        is first carried out on a balanced unit-square multiple of ``alpha``.

    Raises
    ------
    ValueError
        If ``alpha`` is the zero element. The package follows the paper and
        treats exact length only for nonzero inputs.
    """
    if alpha == ZERO:
        raise ValueError("exact_length is defined only for nonzero elements of O_D")
    if not is_representable(order, alpha):
        return None

    p = order.pythagoras_number()
    beta, _ = order.balance_by_pell_unit(alpha) if balance else (alpha, order.one)
    S = relevant_squares(order, beta)
    Hs: Set[Pair] = set(S)

    if beta in Hs:
        return 1

    for s in S:
        diff = order.sub(beta, s)
        if diff in Hs:
            return 2

    if p == 3:
        return 3

    T2 = _pair_sum_set(order, S)

    for s in S:
        diff = order.sub(beta, s)
        if diff in T2:
            return 3

    if p == 4:
        return 4

    for t in T2:
        diff = order.sub(beta, t)
        if diff in T2:
            return 4

    return 5


def decomposition(order: RealQuadraticOrder, alpha: Pair, *, balance: bool = True) -> Optional[List[Pair]]:
    """Return a minimal sum-of-squares decomposition of nonzero ``alpha``.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    alpha : tuple[int, int]
        Nonzero target element in coefficient form.

    Returns
    -------
    list[tuple[int, int]] or None
        A list ``[x_1, ..., x_k]`` of square roots with
        ``alpha = x_1^2 + ... + x_k^2`` and ``k = ell_D(alpha)``, or ``None``
        when ``alpha`` is not representable. When ``balance`` is ``True``
        (the default), decomposition recovery proceeds from a balanced
        unit-square multiple of ``alpha`` and rescales the roots back.

    Raises
    ------
    ValueError
        If ``alpha`` is the zero element. The package follows the paper and
        treats constructive output only for nonzero inputs.
    """
    if alpha == ZERO:
        raise ValueError("decomposition is defined only for nonzero elements of O_D")
    if not is_representable(order, alpha):
        return None

    p = order.pythagoras_number()
    beta, backscale = order.balance_by_pell_unit(alpha) if balance else (alpha, order.one)
    S_list = relevant_squares(order, beta)
    Hs: Set[Pair] = set(S_list)
    sq_to_root = relevant_squares_with_roots(order, beta)

    def _root(square: Pair) -> Pair:
        return sq_to_root[square]

    def _rescale(roots: List[Pair]) -> List[Pair]:
        return [order.mul(backscale, root) for root in roots]

    if beta in Hs:
        return _rescale([_root(beta)])

    for s in S_list:
        diff = order.sub(beta, s)
        if diff in Hs:
            return _rescale([_root(s), _root(diff)])

    T2_witness = _pair_sum_witness(order, S_list)

    for s in S_list:
        diff = order.sub(beta, s)
        if diff in T2_witness:
            w1, w2 = T2_witness[diff]
            return _rescale([_root(s), _root(w1), _root(w2)])

    if p == 3:
        raise RuntimeError("Unreachable: representable element in a ternary case without a 3-square decomposition")

    for t, (w1, w2) in T2_witness.items():
        diff = order.sub(beta, t)
        if diff in T2_witness:
            w3, w4 = T2_witness[diff]
            return _rescale([_root(w1), _root(w2), _root(w3), _root(w4)])

    if p == 4:
        raise RuntimeError("Unreachable: representable element in a quaternary case without a 4-square decomposition")

    for s in sorted(S_list, key=lambda square: (order.trace(square), square[0], square[1]), reverse=True):
        rem = order.sub(beta, s)
        for t, (w1, w2) in T2_witness.items():
            diff = order.sub(rem, t)
            if diff in T2_witness:
                w3, w4 = T2_witness[diff]
                return _rescale([_root(s), _root(w1), _root(w2), _root(w3), _root(w4)])

    raise RuntimeError("Representable element but no 5-square decomposition found")


def _pair_sum_set(order: RealQuadraticOrder, squares: List[Pair]) -> Set[Pair]:
    """Return the unordered pair-sum set ``{s_i + s_j}``.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    squares : list[tuple[int, int]]
        List of relevant squares.

    Returns
    -------
    set[tuple[int, int]]
        Distinct sums of two elements of ``squares``.

    Warns
    -----
    RuntimeWarning
        When the explicit ``T_2`` construction is likely to require substantial
        quadratic memory. In that regime the Schroeppel--Shamir variant
        discussed in the paper may be preferable.
    """
    if len(squares) > 5000:
        warnings.warn(
            "Building T2 explicitly uses O(M^2) memory; for very large search sets consider a lazy meet-in-the-middle variant.",
            RuntimeWarning,
            stacklevel=2,
        )
    result: Set[Pair] = set()
    for index, s1 in enumerate(squares):
        for s2 in squares[index:]:
            result.add(order.add(s1, s2))
    return result


def _pair_sum_witness(order: RealQuadraticOrder, squares: List[Pair]) -> Dict[Pair, Tuple[Pair, Pair]]:
    """Return one witness pair for every unordered pair sum.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    squares : list[tuple[int, int]]
        List of relevant squares.

    Returns
    -------
    dict[tuple[int, int], tuple[tuple[int, int], tuple[int, int]]]
        Dictionary mapping each pair sum to one witness pair.

    Warns
    -----
    RuntimeWarning
        When the explicit witness dictionary is likely to require substantial
        quadratic memory.
    """
    if len(squares) > 5000:
        warnings.warn(
            "Building T2 witnesses explicitly uses O(M^2) memory; for very large search sets consider a lazy meet-in-the-middle variant.",
            RuntimeWarning,
            stacklevel=2,
        )
    witnesses: Dict[Pair, Tuple[Pair, Pair]] = {}
    for index, s1 in enumerate(squares):
        for s2 in squares[index:]:
            t = order.add(s1, s2)
            witnesses.setdefault(t, (s1, s2))
    return witnesses
