r"""Exact diagonal representability over maximal real quadratic orders.

This module is self-contained. It implements exact arithmetic in the maximal
order ``O_D`` of ``Q(sqrt(D))`` and uses it to solve constructive diagonal
representability problems of the form

.. math::

    \alpha = a_1 x_1^2 + \cdots + a_r x_r^2,

where the coefficients ``a_i`` and the target ``alpha`` are totally positive.

The key search primitive is the weighted trace form

.. math::

    q_a(x) = \operatorname{Tr}(a x^2),

which is a positive definite integral binary quadratic form on ``O_D``.
Weighted values are enumerated exactly by row intervals of the trace ellipse.
The implementation uses only integer arithmetic and exact order comparisons via
trace and norm. No floating-point comparisons occur in the core routines.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import product
from statistics import median
from time import perf_counter
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

Pair = Tuple[int, int]


@dataclass(frozen=True)
class RealQuadraticOrder:
    """The maximal order ``O_D`` of ``Q(sqrt(D))``.

    Parameters
    ----------
    D : int
        Squarefree integer at least ``2``.

    Notes
    -----
    Elements are represented by coefficient pairs ``(a, b)`` relative to the
    integral basis ``{1, omega_D}``, where ``omega_D = sqrt(D)`` when
    ``D % 4 in {2, 3}`` and ``omega_D = (1 + sqrt(D)) / 2`` when
    ``D % 4 == 1``.
    """

    D: int

    def __post_init__(self) -> None:
        """Validate the radicand.

        Raises
        ------
        ValueError
            If ``D`` is less than ``2``, a perfect square, or not squarefree.
        """
        if self.D < 2:
            raise ValueError(f"D must be at least 2, got {self.D}")
        s = math.isqrt(self.D)
        if s * s == self.D:
            raise ValueError(f"D = {self.D} is a perfect square")
        if not _is_squarefree(self.D):
            raise ValueError(f"D = {self.D} is not squarefree")

    @property
    def is_mod1(self) -> bool:
        """Return ``True`` when ``D ≡ 1 (mod 4)``."""
        return self.D % 4 == 1

    @property
    def omega_minpoly(self) -> Tuple[int, int]:
        """Return ``(p, q)`` such that ``omega_D^2 = p + q * omega_D``.

        Returns
        -------
        tuple[int, int]
            Coefficients of the multiplication rule for ``omega_D``.
        """
        if self.is_mod1:
            return ((self.D - 1) // 4, 1)
        return (self.D, 0)

    @property
    def field_discriminant(self) -> int:
        """Return the discriminant of ``Q(sqrt(D))``.

        Returns
        -------
        int
            ``D`` if ``D ≡ 1 (mod 4)`` and ``4D`` otherwise.
        """
        return self.D if self.is_mod1 else 4 * self.D

    @property
    def one(self) -> Pair:
        """Return the multiplicative identity.

        Returns
        -------
        tuple[int, int]
            The pair ``(1, 0)``.
        """
        return (1, 0)

    @property
    def zero(self) -> Pair:
        """Return the additive identity.

        Returns
        -------
        tuple[int, int]
            The pair ``(0, 0)``.
        """
        return (0, 0)

    def add(self, x: Pair, y: Pair) -> Pair:
        """Return ``x + y``.

        Parameters
        ----------
        x, y : tuple[int, int]
            Elements of ``O_D`` in coefficient form.

        Returns
        -------
        tuple[int, int]
            The sum ``x + y``.
        """
        return (x[0] + y[0], x[1] + y[1])

    def sub(self, x: Pair, y: Pair) -> Pair:
        """Return ``x - y``.

        Parameters
        ----------
        x, y : tuple[int, int]
            Elements of ``O_D`` in coefficient form.

        Returns
        -------
        tuple[int, int]
            The difference ``x - y``.
        """
        return (x[0] - y[0], x[1] - y[1])

    def neg(self, x: Pair) -> Pair:
        """Return ``-x``.

        Parameters
        ----------
        x : tuple[int, int]
            Element of ``O_D`` in coefficient form.

        Returns
        -------
        tuple[int, int]
            The additive inverse of ``x``.
        """
        return (-x[0], -x[1])

    def mul(self, x: Pair, y: Pair) -> Pair:
        """Return ``x * y``.

        Parameters
        ----------
        x, y : tuple[int, int]
            Elements of ``O_D`` in coefficient form.

        Returns
        -------
        tuple[int, int]
            The product ``x y``.
        """
        a1, b1 = x
        a2, b2 = y
        p, q = self.omega_minpoly
        return (a1 * a2 + b1 * b2 * p, a1 * b2 + a2 * b1 + b1 * b2 * q)

    def sqr(self, x: Pair) -> Pair:
        """Return ``x^2``.

        Parameters
        ----------
        x : tuple[int, int]
            Element of ``O_D`` in coefficient form.

        Returns
        -------
        tuple[int, int]
            The square ``x^2``.
        """
        u, v = x
        p, q = self.omega_minpoly
        return (u * u + v * v * p, 2 * u * v + v * v * q)

    def pow(self, x: Pair, n: int) -> Pair:
        """Return ``x**n`` for ``n >= 0``.

        Parameters
        ----------
        x : tuple[int, int]
            Base element.
        n : int
            Non-negative exponent.

        Returns
        -------
        tuple[int, int]
            The power ``x**n``.

        Raises
        ------
        ValueError
            If ``n`` is negative.
        """
        if n < 0:
            raise ValueError("Exponent must be non-negative")
        result = self.one
        base = x
        exp = n
        while exp:
            if exp & 1:
                result = self.mul(result, base)
            base = self.sqr(base)
            exp >>= 1
        return result

    def conjugate(self, x: Pair) -> Pair:
        """Return the Galois conjugate of ``x``.

        Parameters
        ----------
        x : tuple[int, int]
            Element of ``O_D``.

        Returns
        -------
        tuple[int, int]
            The conjugate of ``x``.
        """
        a, b = x
        if self.is_mod1:
            return (a + b, -b)
        return (a, -b)

    def trace(self, x: Pair) -> int:
        """Return the field trace of ``x``.

        Parameters
        ----------
        x : tuple[int, int]
            Element of ``O_D``.

        Returns
        -------
        int
            The exact trace ``Tr(x)``.
        """
        a, b = x
        return 2 * a + b if self.is_mod1 else 2 * a

    def norm(self, x: Pair) -> int:
        """Return the field norm of ``x``.

        Parameters
        ----------
        x : tuple[int, int]
            Element of ``O_D``.

        Returns
        -------
        int
            The exact norm ``N(x)``.
        """
        a, b = x
        if self.is_mod1:
            return a * a + a * b - ((self.D - 1) // 4) * b * b
        return a * a - self.D * b * b

    def sigma1_sign(self, x: Pair) -> int:
        """Return the sign of the first real embedding of ``x``.

        Parameters
        ----------
        x : tuple[int, int]
            Element of ``O_D``.

        Returns
        -------
        int
            ``-1`` for negative, ``0`` for zero, and ``1`` for positive.
        """
        a, b = x
        if self.is_mod1:
            return _sign_of_p_plus_q_sqrtD(2 * a + b, b, self.D)
        return _sign_of_p_plus_q_sqrtD(a, b, self.D)

    def sigma2_sign(self, x: Pair) -> int:
        """Return the sign of the second real embedding of ``x``.

        Parameters
        ----------
        x : tuple[int, int]
            Element of ``O_D``.

        Returns
        -------
        int
            ``-1`` for negative, ``0`` for zero, and ``1`` for positive.
        """
        a, b = x
        if self.is_mod1:
            return _sign_of_p_plus_q_sqrtD(2 * a + b, -b, self.D)
        return _sign_of_p_plus_q_sqrtD(a, -b, self.D)

    def sigma1_gt(self, x: Pair, y: Pair) -> bool:
        """Return whether ``sigma_1(x) > sigma_1(y)`` exactly.

        Parameters
        ----------
        x, y : tuple[int, int]
            Elements of ``O_D``.

        Returns
        -------
        bool
            Whether the first embedding of ``x - y`` is positive.
        """
        return self.sigma1_sign(self.sub(x, y)) > 0

    def sigma1_ge(self, x: Pair, y: Pair) -> bool:
        """Return whether ``sigma_1(x) >= sigma_1(y)`` exactly.

        Parameters
        ----------
        x, y : tuple[int, int]
            Elements of ``O_D``.

        Returns
        -------
        bool
            Whether the first embedding of ``x - y`` is non-negative.
        """
        return self.sigma1_sign(self.sub(x, y)) >= 0

    def is_nonnegative(self, x: Pair) -> bool:
        """Return whether both real embeddings of ``x`` are non-negative.

        Parameters
        ----------
        x : tuple[int, int]
            Element of ``O_D``.

        Returns
        -------
        bool
            ``True`` if and only if ``Tr(x) >= 0`` and ``N(x) >= 0``.
        """
        return self.trace(x) >= 0 and self.norm(x) >= 0

    def is_positive(self, x: Pair) -> bool:
        """Return whether both real embeddings of ``x`` are strictly positive.

        Parameters
        ----------
        x : tuple[int, int]
            Element of ``O_D``.

        Returns
        -------
        bool
            ``True`` if and only if ``Tr(x) > 0`` and ``N(x) > 0``.
        """
        return self.trace(x) > 0 and self.norm(x) > 0

    def is_totally_positive(self, x: Pair) -> bool:
        """Alias for :meth:`is_positive`.

        Parameters
        ----------
        x : tuple[int, int]
            Element of ``O_D``.

        Returns
        -------
        bool
            ``True`` if and only if ``x`` is totally positive.
        """
        return self.is_positive(x)

    def le(self, x: Pair, y: Pair) -> bool:
        """Return whether ``x <= y`` in both real embeddings.

        Parameters
        ----------
        x, y : tuple[int, int]
            Elements of ``O_D``.

        Returns
        -------
        bool
            ``True`` if and only if ``y - x`` is totally non-negative.
        """
        return self.is_nonnegative(self.sub(y, x))

    def lt(self, x: Pair, y: Pair) -> bool:
        """Return whether ``x < y`` in both real embeddings.

        Parameters
        ----------
        x, y : tuple[int, int]
            Elements of ``O_D``.

        Returns
        -------
        bool
            ``True`` if and only if ``y - x`` is totally positive.
        """
        return self.is_positive(self.sub(y, x))

    def balanced_totally_positive_unit(self) -> Pair:
        """Return a positive norm-one unit derived from Pell's equation.

        Returns
        -------
        tuple[int, int]
            A positive norm-one unit in ``O_D``.
        """
        x, y = _minimal_pell_solution(self.D)
        if self.is_mod1:
            return (x - y, 2 * y)
        return (x, y)

    def balance_by_unit_square(self, x: Pair) -> Tuple[Pair, Pair]:
        """Balance a totally positive element by even powers of a positive unit.

        Parameters
        ----------
        x : tuple[int, int]
            Totally positive element.

        Returns
        -------
        tuple[tuple[int, int], tuple[int, int]]
            A pair ``(beta, scale)`` with ``beta = x * u^(2k)`` and
            ``scale = u^(-k)`` for a suitable integer ``k``.
        """
        if not self.is_totally_positive(x):
            raise ValueError("Only totally positive elements can be balanced")
        unit = self.balanced_totally_positive_unit()
        unit_inv = self.conjugate(unit)
        unit_sq = self.sqr(unit)
        unit_sq_inv = self.sqr(unit_inv)
        beta = x
        k = 0
        while self._ratio_exceeds_square(beta, unit_sq):
            beta = self.mul(beta, unit_sq_inv)
            k -= 1
        while self._inverse_ratio_exceeds_square(beta, unit_sq):
            beta = self.mul(beta, unit_sq)
            k += 1
        if k >= 0:
            scale = self.pow(unit_inv, k)
        else:
            scale = self.pow(unit, -k)
        return beta, scale

    def _ratio_exceeds_square(self, x: Pair, unit_sq: Pair) -> bool:
        """Return whether ``sigma_1(x) / sigma_2(x)`` exceeds ``sigma_1(unit_sq)``.

        Parameters
        ----------
        x : tuple[int, int]
            Totally positive element.
        unit_sq : tuple[int, int]
            Square of a positive norm-one unit.

        Returns
        -------
        bool
            Whether the first-embedding ratio is too large.
        """
        return self.sigma1_gt(x, self.mul(unit_sq, self.conjugate(x)))

    def _inverse_ratio_exceeds_square(self, x: Pair, unit_sq: Pair) -> bool:
        """Return whether ``sigma_2(x) / sigma_1(x)`` exceeds ``sigma_1(unit_sq)``.

        Parameters
        ----------
        x : tuple[int, int]
            Totally positive element.
        unit_sq : tuple[int, int]
            Square of a positive norm-one unit.

        Returns
        -------
        bool
            Whether the inverse first-embedding ratio is too large.
        """
        return self.sigma1_gt(self.conjugate(x), self.mul(unit_sq, x))


@dataclass(frozen=True)
class WeightedTraceForm:
    """The integral binary quadratic form ``q_a(m + n omega_D)``.

    Parameters
    ----------
    A : int
        Coefficient of ``m^2``.
    B : int
        Coefficient of ``m n``.
    C : int
        Coefficient of ``n^2``.
    delta : int
        Positive discriminant parameter ``4 A C - B^2``.
    """

    A: int
    B: int
    C: int
    delta: int

    def value(self, m: int, n: int) -> int:
        """Return ``A m^2 + B m n + C n^2``.

        Parameters
        ----------
        m, n : int
            Integral coordinates.

        Returns
        -------
        int
            The quadratic-form value.
        """
        return self.A * m * m + self.B * m * n + self.C * n * n

    def row_bound(self, trace_bound: int) -> int:
        """Return the maximal absolute row index for ``q_a <= trace_bound``.

        Parameters
        ----------
        trace_bound : int
            Upper bound for the form value.

        Returns
        -------
        int
            The exact upper bound for ``|n|``.
        """
        return math.isqrt((4 * self.A * trace_bound) // self.delta)

    def sharp_box_bounds(self, trace_bound: int) -> Tuple[int, int]:
        """Return the best axis-aligned bounds for the trace ellipse.

        Parameters
        ----------
        trace_bound : int
            Upper bound for the form value.

        Returns
        -------
        tuple[int, int]
            Bounds ``(m_max, n_max)``.
        """
        m_max = math.isqrt((4 * self.C * trace_bound) // self.delta)
        n_max = math.isqrt((4 * self.A * trace_bound) // self.delta)
        return (m_max, n_max)


@dataclass(frozen=True)
class EnumerationStats:
    """Instrumentation for a weighted-search run.

    Parameters
    ----------
    rows_scanned : int
        Number of canonical integral rows examined.
    trace_candidates : int
        Number of canonical non-zero lattice points scanned under the trace
        inequality.
    distinct_values : int
        Number of distinct weighted values accepted after filtering.
    accepted_roots : int
        Total number of accepted non-zero roots, counting both signs.
    elapsed_seconds : float
        Wall-clock runtime.
    """

    rows_scanned: int
    trace_candidates: int
    distinct_values: int
    accepted_roots: int
    elapsed_seconds: float


@dataclass(frozen=True)
class SearchResult:
    """Exact weighted search output.

    Parameters
    ----------
    value_to_root : dict[tuple[int, int], tuple[int, int]]
        One witness root for each distinct weighted value.
    stats : EnumerationStats
        Enumeration instrumentation.
    roots : tuple[tuple[int, int], ...] | None, default=None
        All accepted roots when explicitly requested.
    """

    value_to_root: Dict[Pair, Pair]
    stats: EnumerationStats
    roots: Optional[Tuple[Pair, ...]] = None

    @property
    def values(self) -> List[Pair]:
        """Return the distinct weighted values in sorted order.

        Returns
        -------
        list[tuple[int, int]]
            Sorted list of distinct weighted values.
        """
        result = list(self.value_to_root)
        result.sort()
        return result


@dataclass(frozen=True)
class RepresentationResult:
    """Constructive output of a representability routine.

    Parameters
    ----------
    represented : bool
        Whether the target is represented.
    roots : tuple[tuple[int, int], ...] | None
        One tuple of roots, returned in the original coefficient order when the
        target is represented.
    values : tuple[tuple[int, int], ...] | None
        The corresponding weighted summands, again in the original coefficient
        order.
    preprocessing_seconds : float
        Time spent building the weighted search sets.
    combination_seconds : float
        Time spent in the DP or meet-in-the-middle combination layer.
    elapsed_seconds : float
        End-to-end runtime of the representability routine.
    state_counts : tuple[int, ...], default=()
        Sizes of the DP state sets after each internal layer.
    left_states : int, default=0
        Number of left half-sums for meet-in-the-middle.
    right_states : int, default=0
        Number of right half-sums for meet-in-the-middle.
    distinct_value_counts : tuple[int, ...], default=()
        Sizes ``|S_{a_i}(alpha)|`` in the internal coefficient order actually
        used by the combination phase.
    coefficient_order : tuple[int, ...], default=()
        Permutation describing the internal coefficient order. Entry ``j`` is
        the original index of the ``j``th internal coefficient.
    left_prefix_counts, right_prefix_counts : tuple[int, ...], default=()
        Prefix half-sum counts for the two meet-in-the-middle halves.
    unique_preprocessed_coeffs : int, default=0
        Number of distinct coefficients for which weighted-value sets were
        actually computed.
    """

    represented: bool
    roots: Optional[Tuple[Pair, ...]]
    values: Optional[Tuple[Pair, ...]]
    preprocessing_seconds: float
    combination_seconds: float
    elapsed_seconds: float
    state_counts: Tuple[int, ...] = ()
    left_states: int = 0
    right_states: int = 0
    distinct_value_counts: Tuple[int, ...] = ()
    coefficient_order: Tuple[int, ...] = ()
    left_prefix_counts: Tuple[int, ...] = ()
    right_prefix_counts: Tuple[int, ...] = ()
    unique_preprocessed_coeffs: int = 0


@dataclass(frozen=True)
class BoundedRepresentabilityResult:
    """Batched bounded representability output.

    Parameters
    ----------
    representables : tuple[tuple[int, int], ...]
        All non-zero represented targets of trace at most ``B``.
    state_counts : tuple[int, ...]
        Sizes of the bounded state sets after each internal coefficient.
    preprocessing_seconds : float
        Time spent building the bounded weighted-value sets.
    combination_seconds : float
        Time spent in the batched dynamic-programming layer.
    elapsed_seconds : float
        End-to-end runtime of the bounded representability routine.
    value_counts : tuple[int, ...]
        Sizes ``|W_i(B)| - 1`` of the non-zero weighted-value sets in the
        internal coefficient order actually used by the combination phase.
    coefficient_order : tuple[int, ...]
        Permutation describing the internal coefficient order.
    unique_preprocessed_coeffs : int
        Number of distinct coefficients for which bounded weighted-value sets
        were actually computed.
    """

    representables: Tuple[Pair, ...]
    state_counts: Tuple[int, ...]
    preprocessing_seconds: float
    combination_seconds: float
    elapsed_seconds: float
    value_counts: Tuple[int, ...]
    coefficient_order: Tuple[int, ...]
    unique_preprocessed_coeffs: int



@dataclass(frozen=True)
class PreparedSingleTarget:
    """Precomputed single-target weighted-value data."""

    witness_maps: Tuple[Dict[Pair, Pair], ...]
    value_lists: Tuple[Tuple[Pair, ...], ...]
    distinct_value_counts: Tuple[int, ...]
    coefficient_order: Tuple[int, ...]
    preprocessing_seconds: float
    unique_preprocessed_coeffs: int


@dataclass(frozen=True)
class PreparedBoundedValues:
    """Precomputed bounded weighted-value data."""

    value_lists: Tuple[Tuple[Pair, ...], ...]
    value_counts: Tuple[int, ...]
    coefficient_order: Tuple[int, ...]
    preprocessing_seconds: float
    unique_preprocessed_coeffs: int


@dataclass(frozen=True)
class EnumerationBenchmarkRow:
    """One benchmark row for weighted enumeration."""

    D: int
    coeff: Pair
    alpha: Pair
    distinct_values: int
    accepted_roots: int
    sharp_box_candidates: int
    crude_box_candidates: int
    row_candidates: int
    sharp_box_ms: float
    crude_box_ms: float
    row_ms: float


@dataclass(frozen=True)
class RepresentabilityBenchmarkRow:
    """One benchmark row for representability routines."""

    D: int
    coeffs: Tuple[Pair, ...]
    alpha: Pair
    coefficient_order: Tuple[int, ...]
    distinct_value_counts: Tuple[int, ...]
    dp_state_counts: Tuple[int, ...]
    mitm_left_states: int
    mitm_right_states: int
    preprocessing_ms: float
    dp_combination_ms: float
    dp_total_ms: float
    mitm_combination_ms: float
    mitm_total_ms: float


@dataclass(frozen=True)
class GenericBaselineBenchmarkRow:
    """One benchmark row for the generic exact representability baseline."""

    D: int
    coeffs: Tuple[Pair, ...]
    alpha: Pair
    coefficient_order: Tuple[int, ...]
    distinct_value_counts: Tuple[int, ...]
    cartesian_product_size: int
    preprocessing_ms: float
    dp_total_ms: float
    mitm_total_ms: float
    generic_total_ms: float


@dataclass(frozen=True)
class BoundedBenchmarkRow:
    """One benchmark row for bounded truant routines."""

    D: int
    coeffs: Tuple[Pair, ...]
    trace_bound: int
    representables: int
    truants: int
    batched_state_counts: Tuple[int, ...]
    value_counts: Tuple[int, ...]
    batched_preprocessing_ms: float
    batched_combination_ms: float
    batched_total_ms: float
    naive_total_ms: float


@dataclass(frozen=True)
class OptimisationBenchmarkRow:
    """One benchmark row for coefficient-order optimisation and caching."""

    D: int
    coeffs: Tuple[Pair, ...]
    alpha: Pair
    original_distinct_value_counts: Tuple[int, ...]
    optimised_distinct_value_counts: Tuple[int, ...]
    original_state_counts: Tuple[int, ...]
    optimised_state_counts: Tuple[int, ...]
    coefficient_order: Tuple[int, ...]
    uncached_preprocessing_ms: float
    cached_preprocessing_ms: float
    original_total_ms: float
    optimised_total_ms: float


@dataclass(frozen=True)
class BalancingBenchmarkRow:
    """One deterministic row illustrating unit-square balancing."""

    D: int
    coeff: Pair
    balanced_coeff: Pair
    trace_bound: int
    original_row_candidates: int
    balanced_row_candidates: int
    original_sharp_box_candidates: int
    balanced_sharp_box_candidates: int


def weighted_trace_form(order: RealQuadraticOrder, coeff: Pair) -> WeightedTraceForm:
    """Return the weighted trace form attached to ``coeff``.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    coeff : tuple[int, int]
        Totally positive coefficient ``a``.

    Returns
    -------
    WeightedTraceForm
        The integral form ``q_a(m + n omega_D)``.

    Raises
    ------
    ValueError
        If ``coeff`` is not totally positive.
    """
    if not order.is_totally_positive(coeff):
        raise ValueError("Weighted coefficient must be totally positive")
    e1 = (1, 0)
    e2 = (0, 1)
    A = order.trace(order.mul(coeff, order.sqr(e1)))
    C = order.trace(order.mul(coeff, order.sqr(e2)))
    e12 = order.add(e1, e2)
    q12 = order.trace(order.mul(coeff, order.sqr(e12)))
    B = q12 - A - C
    delta = 4 * A * C - B * B
    if A <= 0 or C <= 0 or delta <= 0:
        raise ValueError("The weighted trace form must be positive definite")
    return WeightedTraceForm(A=A, B=B, C=C, delta=delta)


def exact_row_interval(form: WeightedTraceForm, trace_bound: int, n: int) -> Optional[Tuple[int, int]]:
    """Return the exact integer interval on row ``n``.

    Parameters
    ----------
    form : WeightedTraceForm
        Positive definite form ``A m^2 + B m n + C n^2``.
    trace_bound : int
        Bound ``T`` in the inequality ``q_a(m, n) <= T``.
    n : int
        Fixed row index.

    Returns
    -------
    tuple[int, int] | None
        The exact interval ``[m_min, m_max]`` or ``None`` if the row is empty.
    """
    disc = 4 * form.A * trace_bound - form.delta * n * n
    if disc < 0:
        return None
    s = math.isqrt(disc)
    left = _ceil_div(-form.B * n - s, 2 * form.A)
    right = _floor_div(-form.B * n + s, 2 * form.A)
    if left > right:
        return None
    return (left, right)


def discriminant_identity(order: RealQuadraticOrder, coeff: Pair) -> int:
    """Return ``4 * Disc(K_D) * N(coeff)``.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    coeff : tuple[int, int]
        Totally positive coefficient.

    Returns
    -------
    int
        The predicted discriminant parameter of ``q_a``.
    """
    return 4 * order.field_discriminant * order.norm(coeff)


def _canonical_m_values(m_min: int, m_max: int, n: int) -> range:
    """Return the canonical integral abscissae on a fixed row.

    Parameters
    ----------
    m_min, m_max : int
        Endpoints of the exact row interval.
    n : int
        Row index.

    Returns
    -------
    range
        The canonical abscissae. For positive rows all abscissae are retained;
        on the central row only positive abscissae are kept.
    """
    if n == 0:
        if m_max < 1:
            return range(0, 0)
        return range(max(1, m_min), m_max + 1)
    return range(m_min, m_max + 1)


def _append_all_sign_roots(roots: List[Pair], x: Pair) -> None:
    """Append both sign choices of a non-zero root.

    Parameters
    ----------
    roots : list[tuple[int, int]]
        Output list.
    x : tuple[int, int]
        Canonical non-zero root.
    """
    roots.append(x)
    roots.append((-x[0], -x[1]))


def enumerate_weighted_values_by_trace(
    order: RealQuadraticOrder,
    coeff: Pair,
    trace_bound: int,
    return_all_roots: bool = False,
) -> SearchResult:
    """Enumerate all non-zero weighted values of trace at most ``trace_bound``.

    The routine scans only a canonical half of the trace ellipse, namely the
    roots with ``n > 0`` together with the roots on the central row with
    ``m > 0``. This visits exactly one representative of each non-zero
    ``±``-pair.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    coeff : tuple[int, int]
        Totally positive coefficient.
    trace_bound : int
        Upper trace bound.
    return_all_roots : bool, default=False
        Whether to store all accepted roots explicitly.

    Returns
    -------
    SearchResult
        Distinct values, witness roots, and instrumentation.
    """
    if trace_bound < 0:
        raise ValueError("Trace bound must be non-negative")
    form = weighted_trace_form(order, coeff)
    n_max = form.row_bound(trace_bound)
    value_to_root: Dict[Pair, Pair] = {}
    roots: Optional[List[Pair]] = [] if return_all_roots else None
    rows_scanned = 0
    trace_candidates = 0

    start = perf_counter()
    for n in range(0, n_max + 1):
        rows_scanned += 1
        interval = exact_row_interval(form, trace_bound, n)
        if interval is None:
            continue
        m_min, m_max = interval
        canonical_ms = _canonical_m_values(m_min, m_max, n)
        trace_candidates += len(canonical_ms)
        for m in canonical_ms:
            x = (m, n)
            value = order.mul(coeff, order.sqr(x))
            value_to_root.setdefault(value, x)
            if roots is not None:
                _append_all_sign_roots(roots, x)
    elapsed = perf_counter() - start

    distinct_values = len(value_to_root)
    accepted_roots = len(roots) if roots is not None else 2 * distinct_values
    roots_tuple: Optional[Tuple[Pair, ...]]
    if roots is None:
        roots_tuple = None
    else:
        roots.sort()
        roots_tuple = tuple(roots)
    stats = EnumerationStats(
        rows_scanned=rows_scanned,
        trace_candidates=trace_candidates,
        distinct_values=distinct_values,
        accepted_roots=accepted_roots,
        elapsed_seconds=elapsed,
    )
    return SearchResult(value_to_root=value_to_root, stats=stats, roots=roots_tuple)


def enumerate_weighted_search(
    order: RealQuadraticOrder,
    coeff: Pair,
    alpha: Pair,
    return_all_roots: bool = False,
) -> SearchResult:
    """Enumerate the weighted search set ``S_a(alpha)`` exactly.

    The routine scans only a canonical half of the trace ellipse and keeps one
    witness root for each distinct weighted value.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    coeff : tuple[int, int]
        Totally positive coefficient.
    alpha : tuple[int, int]
        Totally positive target.
    return_all_roots : bool, default=False
        Whether to store all accepted roots explicitly.

    Returns
    -------
    SearchResult
        Distinct weighted values ``a x^2`` with ``0 < a x^2 <= alpha``.
    """
    if not order.is_totally_positive(alpha):
        raise ValueError("Target must be totally positive")
    form = weighted_trace_form(order, coeff)
    trace_bound = order.trace(alpha)
    n_max = form.row_bound(trace_bound)
    value_to_root: Dict[Pair, Pair] = {}
    roots: Optional[List[Pair]] = [] if return_all_roots else None
    rows_scanned = 0
    trace_candidates = 0

    start = perf_counter()
    for n in range(0, n_max + 1):
        rows_scanned += 1
        interval = exact_row_interval(form, trace_bound, n)
        if interval is None:
            continue
        m_min, m_max = interval
        canonical_ms = _canonical_m_values(m_min, m_max, n)
        trace_candidates += len(canonical_ms)
        for m in canonical_ms:
            x = (m, n)
            value = order.mul(coeff, order.sqr(x))
            if order.le(value, alpha):
                value_to_root.setdefault(value, x)
                if roots is not None:
                    _append_all_sign_roots(roots, x)
    elapsed = perf_counter() - start

    distinct_values = len(value_to_root)
    accepted_roots = len(roots) if roots is not None else 2 * distinct_values
    roots_tuple: Optional[Tuple[Pair, ...]]
    if roots is None:
        roots_tuple = None
    else:
        roots.sort()
        roots_tuple = tuple(roots)
    stats = EnumerationStats(
        rows_scanned=rows_scanned,
        trace_candidates=trace_candidates,
        distinct_values=distinct_values,
        accepted_roots=accepted_roots,
        elapsed_seconds=elapsed,
    )
    return SearchResult(value_to_root=value_to_root, stats=stats, roots=roots_tuple)


def enumerate_weighted_search_box(
    order: RealQuadraticOrder,
    coeff: Pair,
    alpha: Pair,
    use_sharp_box: bool = True,
    return_all_roots: bool = False,
) -> SearchResult:
    """Enumerate the weighted search set by canonical box-and-filter.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    coeff : tuple[int, int]
        Totally positive coefficient.
    alpha : tuple[int, int]
        Totally positive target.
    use_sharp_box : bool, default=True
        Whether to use the sharp box from completing the square. If ``False``,
        use the deliberately crude box.
    return_all_roots : bool, default=False
        Whether to store all accepted roots explicitly.

    Returns
    -------
    SearchResult
        Search result and instrumentation.
    """
    if not order.is_totally_positive(coeff):
        raise ValueError("Weighted coefficient must be totally positive")
    if not order.is_totally_positive(alpha):
        raise ValueError("Target must be totally positive")
    form = weighted_trace_form(order, coeff)
    trace_bound = order.trace(alpha)
    if use_sharp_box:
        m_max, n_max = form.sharp_box_bounds(trace_bound)
    else:
        m_max = math.isqrt(4 * form.C * trace_bound)
        n_max = math.isqrt(4 * form.A * trace_bound)

    value_to_root: Dict[Pair, Pair] = {}
    roots: Optional[List[Pair]] = [] if return_all_roots else None
    rows_scanned = n_max + 1
    trace_candidates = 0

    start = perf_counter()
    for n in range(0, n_max + 1):
        m_values = range(1, m_max + 1) if n == 0 else range(-m_max, m_max + 1)
        trace_candidates += len(m_values)
        for m in m_values:
            if form.value(m, n) > trace_bound:
                continue
            x = (m, n)
            value = order.mul(coeff, order.sqr(x))
            if order.le(value, alpha):
                value_to_root.setdefault(value, x)
                if roots is not None:
                    _append_all_sign_roots(roots, x)
    elapsed = perf_counter() - start

    distinct_values = len(value_to_root)
    accepted_roots = len(roots) if roots is not None else 2 * distinct_values
    roots_tuple: Optional[Tuple[Pair, ...]]
    if roots is None:
        roots_tuple = None
    else:
        roots.sort()
        roots_tuple = tuple(roots)
    stats = EnumerationStats(
        rows_scanned=rows_scanned,
        trace_candidates=trace_candidates,
        distinct_values=distinct_values,
        accepted_roots=accepted_roots,
        elapsed_seconds=elapsed,
    )
    return SearchResult(value_to_root=value_to_root, stats=stats, roots=roots_tuple)


def _optimised_coefficient_order(
    coeffs: Sequence[Pair],
    counts: Sequence[int],
    reorder_by_value_count: bool,
) -> Tuple[int, ...]:
    """Return the internal coefficient order used by the algorithms.

    Parameters
    ----------
    coeffs : sequence[tuple[int, int]]
        Original diagonal coefficients.
    counts : sequence[int]
        Search-set sizes associated with the coefficients.
    reorder_by_value_count : bool
        Whether to sort by non-decreasing search-set size.

    Returns
    -------
    tuple[int, ...]
        The chosen internal order. Entry ``j`` is the original index of the
        ``j``th internal coefficient.
    """
    if not reorder_by_value_count:
        return tuple(range(len(coeffs)))
    return tuple(sorted(range(len(coeffs)), key=lambda i: (counts[i], coeffs[i], i)))


def _prepare_single_target(
    order: RealQuadraticOrder,
    coeffs: Sequence[Pair],
    alpha: Pair,
    reorder_by_value_count: bool = True,
    cache_repeated_coeffs: bool = True,
) -> PreparedSingleTarget:
    """Precompute the weighted-value sets for one target.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    coeffs : sequence[tuple[int, int]]
        Totally positive diagonal coefficients.
    alpha : tuple[int, int]
        Totally positive target.
    reorder_by_value_count : bool, default=True
        Whether to sort coefficients by increasing ``|S_{a_i}(alpha)|`` before
        the combination phase.
    cache_repeated_coeffs : bool, default=True
        Whether to preprocess each distinct repeated coefficient only once.

    Returns
    -------
    PreparedSingleTarget
        Witness maps, value lists including zero, internal coefficient order,
        and preprocessing instrumentation.
    """
    zero = order.zero
    witness_maps_raw: List[Dict[Pair, Pair]] = []
    value_lists_raw: List[Tuple[Pair, ...]] = []
    distinct_value_counts_raw: List[int] = []
    cache: Dict[Pair, SearchResult] = {}
    unique_preprocessed_coeffs = 0
    start = perf_counter()
    for coeff in coeffs:
        if cache_repeated_coeffs and coeff in cache:
            search = cache[coeff]
        else:
            search = enumerate_weighted_search(order, coeff, alpha, return_all_roots=False)
            if cache_repeated_coeffs:
                cache[coeff] = search
            unique_preprocessed_coeffs += 1
        witness_maps_raw.append(search.value_to_root)
        value_lists_raw.append(tuple([zero] + search.values))
        distinct_value_counts_raw.append(search.stats.distinct_values)
    elapsed = perf_counter() - start
    if not cache_repeated_coeffs:
        unique_preprocessed_coeffs = len(coeffs)
    coefficient_order = _optimised_coefficient_order(
        coeffs,
        distinct_value_counts_raw,
        reorder_by_value_count=reorder_by_value_count,
    )
    witness_maps = tuple(witness_maps_raw[i] for i in coefficient_order)
    value_lists = tuple(value_lists_raw[i] for i in coefficient_order)
    distinct_value_counts = tuple(distinct_value_counts_raw[i] for i in coefficient_order)
    return PreparedSingleTarget(
        witness_maps=witness_maps,
        value_lists=value_lists,
        distinct_value_counts=distinct_value_counts,
        coefficient_order=coefficient_order,
        preprocessing_seconds=elapsed,
        unique_preprocessed_coeffs=unique_preprocessed_coeffs,
    )


def _dp_from_value_lists(
    order: RealQuadraticOrder,
    value_lists: Sequence[Sequence[Pair]],
    alpha: Pair,
) -> Tuple[bool, Optional[Tuple[Pair, ...]], Tuple[int, ...], float]:
    """Run the constructive DP layer from precomputed value lists.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    value_lists : sequence[sequence[tuple[int, int]]]
        Candidate value lists for each coefficient.
    alpha : tuple[int, int]
        Totally positive target.

    Returns
    -------
    tuple
        Representability flag, chosen values when represented, state counts,
        and the combination time.
    """
    zero = order.zero
    partial: Dict[Pair, None] = {zero: None}
    predecessors: List[Dict[Pair, Tuple[Pair, Pair]]] = []
    state_counts: List[int] = []
    start = perf_counter()
    for values in value_lists:
        next_partial: Dict[Pair, None] = {}
        pred: Dict[Pair, Tuple[Pair, Pair]] = {}
        for y in partial:
            for s in values:
                t = order.add(y, s)
                if order.le(t, alpha):
                    if t not in next_partial:
                        next_partial[t] = None
                        pred[t] = (y, s)
        partial = next_partial
        predecessors.append(pred)
        state_counts.append(len(partial))
    elapsed = perf_counter() - start

    if alpha not in partial:
        return (False, None, tuple(state_counts), elapsed)

    chosen_values: List[Pair] = []
    current = alpha
    for level in range(len(value_lists) - 1, -1, -1):
        prev, summand = predecessors[level][current]
        chosen_values.append(summand)
        current = prev
    chosen_values.reverse()
    return (True, tuple(chosen_values), tuple(state_counts), elapsed)

def diagonal_representability_dp(
    order: RealQuadraticOrder,
    coeffs: Sequence[Pair],
    alpha: Pair,
    reorder_by_value_count: bool = True,
    cache_repeated_coeffs: bool = True,
) -> RepresentationResult:
    """Decide diagonal representability by dynamic programming.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    coeffs : sequence[tuple[int, int]]
        Totally positive diagonal coefficients.
    alpha : tuple[int, int]
        Totally positive target.
    reorder_by_value_count : bool, default=True
        Whether to sort coefficients by increasing search-set size before the
        combination phase.
    cache_repeated_coeffs : bool, default=True
        Whether to preprocess each distinct repeated coefficient only once.

    Returns
    -------
    RepresentationResult
        Exact decision, witness roots in the original coefficient order, state
        counts, and split timings.
    """
    if not order.is_totally_positive(alpha):
        raise ValueError("Target must be totally positive")
    zero = order.zero
    prepared = _prepare_single_target(
        order,
        coeffs,
        alpha,
        reorder_by_value_count=reorder_by_value_count,
        cache_repeated_coeffs=cache_repeated_coeffs,
    )
    represented, chosen_values_internal, state_counts, combination_seconds = _dp_from_value_lists(
        order,
        prepared.value_lists,
        alpha,
    )
    elapsed = prepared.preprocessing_seconds + combination_seconds

    if not represented or chosen_values_internal is None:
        return RepresentationResult(
            represented=False,
            roots=None,
            values=None,
            preprocessing_seconds=prepared.preprocessing_seconds,
            combination_seconds=combination_seconds,
            elapsed_seconds=elapsed,
            state_counts=state_counts,
            distinct_value_counts=prepared.distinct_value_counts,
            coefficient_order=prepared.coefficient_order,
            unique_preprocessed_coeffs=prepared.unique_preprocessed_coeffs,
        )

    chosen_values_original: List[Pair] = [zero for _ in coeffs]
    roots_original: List[Pair] = [zero for _ in coeffs]
    for internal_index, original_index in enumerate(prepared.coefficient_order):
        summand = chosen_values_internal[internal_index]
        chosen_values_original[original_index] = summand
        if summand == zero:
            roots_original[original_index] = zero
        else:
            roots_original[original_index] = prepared.witness_maps[internal_index][summand]
    return RepresentationResult(
        represented=True,
        roots=tuple(roots_original),
        values=tuple(chosen_values_original),
        preprocessing_seconds=prepared.preprocessing_seconds,
        combination_seconds=combination_seconds,
        elapsed_seconds=elapsed,
        state_counts=state_counts,
        distinct_value_counts=prepared.distinct_value_counts,
        coefficient_order=prepared.coefficient_order,
        unique_preprocessed_coeffs=prepared.unique_preprocessed_coeffs,
    )


def diagonal_representability_mitm(
    order: RealQuadraticOrder,
    coeffs: Sequence[Pair],
    alpha: Pair,
    reorder_by_value_count: bool = True,
    cache_repeated_coeffs: bool = True,
) -> RepresentationResult:
    """Decide diagonal representability by meet-in-the-middle.

    The implementation splits after the first ``floor(r / 2)`` coefficients in
    the internally reordered sequence.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    coeffs : sequence[tuple[int, int]]
        Totally positive diagonal coefficients.
    alpha : tuple[int, int]
        Totally positive target.
    reorder_by_value_count : bool, default=True
        Whether to sort coefficients by increasing search-set size before the
        combination phase.
    cache_repeated_coeffs : bool, default=True
        Whether to preprocess each distinct repeated coefficient only once.

    Returns
    -------
    RepresentationResult
        Exact decision, witness roots in the original coefficient order, half-
        sum sizes, and split timings.
    """
    if not order.is_totally_positive(alpha):
        raise ValueError("Target must be totally positive")
    zero = order.zero
    prepared = _prepare_single_target(
        order,
        coeffs,
        alpha,
        reorder_by_value_count=reorder_by_value_count,
        cache_repeated_coeffs=cache_repeated_coeffs,
    )

    split = len(coeffs) // 2
    left_values = prepared.value_lists[:split]
    right_values = prepared.value_lists[split:]

    start = perf_counter()
    left_sums, left_prefix_counts = _half_sums(order, left_values, alpha)
    right_sums, right_prefix_counts = _half_sums(order, right_values, alpha)
    if len(left_sums) <= len(right_sums):
        scan_left = True
        small_sums = left_sums
        large_sums = right_sums
    else:
        scan_left = False
        small_sums = right_sums
        large_sums = left_sums

    chosen_values_internal: Optional[Tuple[Pair, ...]] = None
    for partial_sum, partial_choice in small_sums.items():
        need = order.sub(alpha, partial_sum)
        other_choice = large_sums.get(need)
        if other_choice is None:
            continue
        if scan_left:
            chosen_values_internal = tuple(partial_choice) + tuple(other_choice)
        else:
            chosen_values_internal = tuple(other_choice) + tuple(partial_choice)
        break
    combination_seconds = perf_counter() - start
    elapsed = prepared.preprocessing_seconds + combination_seconds

    if chosen_values_internal is None:
        return RepresentationResult(
            represented=False,
            roots=None,
            values=None,
            preprocessing_seconds=prepared.preprocessing_seconds,
            combination_seconds=combination_seconds,
            elapsed_seconds=elapsed,
            left_states=len(left_sums),
            right_states=len(right_sums),
            distinct_value_counts=prepared.distinct_value_counts,
            coefficient_order=prepared.coefficient_order,
            left_prefix_counts=left_prefix_counts,
            right_prefix_counts=right_prefix_counts,
            unique_preprocessed_coeffs=prepared.unique_preprocessed_coeffs,
        )

    chosen_values_original: List[Pair] = [zero for _ in coeffs]
    roots_original: List[Pair] = [zero for _ in coeffs]
    for internal_index, original_index in enumerate(prepared.coefficient_order):
        summand = chosen_values_internal[internal_index]
        chosen_values_original[original_index] = summand
        if summand == zero:
            roots_original[original_index] = zero
        else:
            roots_original[original_index] = prepared.witness_maps[internal_index][summand]
    return RepresentationResult(
        represented=True,
        roots=tuple(roots_original),
        values=tuple(chosen_values_original),
        preprocessing_seconds=prepared.preprocessing_seconds,
        combination_seconds=combination_seconds,
        elapsed_seconds=elapsed,
        left_states=len(left_sums),
        right_states=len(right_sums),
        distinct_value_counts=prepared.distinct_value_counts,
        coefficient_order=prepared.coefficient_order,
        left_prefix_counts=left_prefix_counts,
        right_prefix_counts=right_prefix_counts,
        unique_preprocessed_coeffs=prepared.unique_preprocessed_coeffs,
    )


def _half_sums(
    order: RealQuadraticOrder,
    value_lists: Sequence[Sequence[Pair]],
    alpha: Pair,
) -> Tuple[Dict[Pair, Tuple[Pair, ...]], Tuple[int, ...]]:
    """Return one choice tuple for each bounded half-sum.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    value_lists : sequence[sequence[tuple[int, int]]]
        Candidate values for each coefficient in the half-lattice.
    alpha : tuple[int, int]
        Bounding target.

    Returns
    -------
    tuple
        A dictionary storing one choice tuple for each reachable half-sum, and
        the prefix half-sum counts after each layer.
    """
    states: Dict[Pair, Tuple[Pair, ...]] = {order.zero: tuple()}
    prefix_counts: List[int] = []
    for values in value_lists:
        next_states: Dict[Pair, Tuple[Pair, ...]] = {}
        for partial_sum, choice in states.items():
            for summand in values:
                total = order.add(partial_sum, summand)
                if order.le(total, alpha):
                    next_states.setdefault(total, choice + (summand,))
        states = next_states
        prefix_counts.append(len(states))
    return states, tuple(prefix_counts)


def _cartesian_product_from_value_lists(
    order: RealQuadraticOrder,
    value_lists: Sequence[Sequence[Pair]],
    alpha: Pair,
) -> Tuple[bool, Optional[Tuple[Pair, ...]], int, float]:
    """Solve representability by explicit Cartesian-product search.

    This routine is intended as a reproducible generic exact baseline for the
    specialised DP and meet-in-the-middle layers. It uses the same precomputed
    weighted value lists but expands the full Cartesian product instead of
    exploiting partial-sum reuse.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    value_lists : sequence[sequence[tuple[int, int]]]
        Candidate value lists, usually including zero.
    alpha : tuple[int, int]
        Totally positive target.

    Returns
    -------
    tuple
        Representability flag, one chosen tuple of values when represented, the
        total number of Cartesian-product combinations visited, and the elapsed
        combination time.
    """
    cartesian_product_size = 1
    for values in value_lists:
        cartesian_product_size *= len(values)

    chosen_values: Optional[Tuple[Pair, ...]] = None
    start = perf_counter()
    for choice in product(*value_lists):
        total = order.zero
        for summand in choice:
            total = order.add(total, summand)
        if total == alpha and chosen_values is None:
            chosen_values = tuple(choice)
    elapsed = perf_counter() - start
    return (chosen_values is not None, chosen_values, cartesian_product_size, elapsed)


def enumerate_totally_positive_targets(order: RealQuadraticOrder, trace_bound: int) -> List[Pair]:
    """Enumerate all non-zero totally positive elements of bounded trace.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    trace_bound : int
        Upper bound for the trace.

    Returns
    -------
    list[tuple[int, int]]
        Sorted list of all non-zero totally positive elements ``alpha`` with
        ``Tr(alpha) <= trace_bound``.
    """
    if trace_bound < 0:
        raise ValueError("Trace bound must be non-negative")
    result: List[Pair] = []
    if order.is_mod1:
        for n in range(-trace_bound, trace_bound + 1):
            lower_trace = math.isqrt(order.D * n * n) + 1
            t0 = lower_trace if (lower_trace - n) % 2 == 0 else lower_trace + 1
            for tr in range(t0, trace_bound + 1, 2):
                m = (tr - n) // 2
                alpha = (m, n)
                if order.is_totally_positive(alpha):
                    result.append(alpha)
    else:
        m_max = trace_bound // 2
        for n in range(-m_max, m_max + 1):
            lower_m = math.isqrt(order.D * n * n) + 1
            for m in range(lower_m, m_max + 1):
                alpha = (m, n)
                if order.is_totally_positive(alpha):
                    result.append(alpha)
    result.sort(key=lambda x: (order.trace(x), x[0], x[1]))
    return result


def _prepare_bounded_value_lists(
    order: RealQuadraticOrder,
    coeffs: Sequence[Pair],
    trace_bound: int,
    reorder_by_value_count: bool = True,
    cache_repeated_coeffs: bool = True,
) -> PreparedBoundedValues:
    """Precompute the bounded weighted-value sets.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    coeffs : sequence[tuple[int, int]]
        Totally positive diagonal coefficients.
    trace_bound : int
        Upper trace bound.
    reorder_by_value_count : bool, default=True
        Whether to sort coefficients by increasing ``|W_i(B)|`` before the
        combination phase.
    cache_repeated_coeffs : bool, default=True
        Whether to preprocess each distinct repeated coefficient only once.

    Returns
    -------
    PreparedBoundedValues
        Bounded value lists including zero, internal coefficient order, and
        preprocessing instrumentation.
    """
    zero = order.zero
    value_lists_raw: List[Tuple[Pair, ...]] = []
    value_counts_raw: List[int] = []
    cache: Dict[Pair, SearchResult] = {}
    unique_preprocessed_coeffs = 0
    start = perf_counter()
    for coeff in coeffs:
        if cache_repeated_coeffs and coeff in cache:
            search = cache[coeff]
        else:
            search = enumerate_weighted_values_by_trace(order, coeff, trace_bound, return_all_roots=False)
            if cache_repeated_coeffs:
                cache[coeff] = search
            unique_preprocessed_coeffs += 1
        value_lists_raw.append(tuple([zero] + search.values))
        value_counts_raw.append(search.stats.distinct_values)
    elapsed = perf_counter() - start
    if not cache_repeated_coeffs:
        unique_preprocessed_coeffs = len(coeffs)
    coefficient_order = _optimised_coefficient_order(
        coeffs,
        value_counts_raw,
        reorder_by_value_count=reorder_by_value_count,
    )
    value_lists = tuple(value_lists_raw[i] for i in coefficient_order)
    value_counts = tuple(value_counts_raw[i] for i in coefficient_order)
    return PreparedBoundedValues(
        value_lists=value_lists,
        value_counts=value_counts,
        coefficient_order=coefficient_order,
        preprocessing_seconds=elapsed,
        unique_preprocessed_coeffs=unique_preprocessed_coeffs,
    )


def _batched_bounded_from_value_lists(
    order: RealQuadraticOrder,
    value_lists: Sequence[Sequence[Pair]],
    trace_bound: int,
) -> Tuple[Tuple[Pair, ...], Tuple[int, ...], float]:
    """Run the bounded batched dynamic program from precomputed value lists.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    value_lists : sequence[sequence[tuple[int, int]]]
        Bounded weighted-value lists including zero.
    trace_bound : int
        Upper trace bound.

    Returns
    -------
    tuple
        Represented non-zero targets, state counts, and combination time.

    Notes
    -----
    The distinguished zero state is retained at every intermediate stage. It
    represents the choice in which all coefficients processed so far contribute
    the zero summand.
    """
    zero = order.zero
    current = {zero}
    state_counts: List[int] = []
    start = perf_counter()
    for values in value_lists:
        next_states = set()
        for y in current:
            for s in values:
                t = order.add(y, s)
                if t == zero or (order.trace(t) <= trace_bound and order.is_totally_positive(t)):
                    next_states.add(t)
        current = next_states
        state_counts.append(len(current))
    elapsed = perf_counter() - start
    representables = sorted((x for x in current if x != zero), key=lambda z: (order.trace(z), z[0], z[1]))
    return tuple(representables), tuple(state_counts), elapsed


def batched_bounded_representables(
    order: RealQuadraticOrder,
    coeffs: Sequence[Pair],
    trace_bound: int,
    reorder_by_value_count: bool = True,
    cache_repeated_coeffs: bool = True,
) -> BoundedRepresentabilityResult:
    """Compute all bounded represented targets at once.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    coeffs : sequence[tuple[int, int]]
        Totally positive diagonal coefficients.
    trace_bound : int
        Upper trace bound ``B``.
    reorder_by_value_count : bool, default=True
        Whether to sort coefficients by increasing ``|W_i(B)|`` before the
        combination phase.
    cache_repeated_coeffs : bool, default=True
        Whether to preprocess each distinct repeated coefficient only once.

    Returns
    -------
    BoundedRepresentabilityResult
        All represented non-zero targets of trace at most ``B``.
    """
    if trace_bound < 0:
        raise ValueError("Trace bound must be non-negative")
    prepared = _prepare_bounded_value_lists(
        order,
        coeffs,
        trace_bound,
        reorder_by_value_count=reorder_by_value_count,
        cache_repeated_coeffs=cache_repeated_coeffs,
    )
    representables, state_counts, combination_seconds = _batched_bounded_from_value_lists(
        order,
        prepared.value_lists,
        trace_bound,
    )
    elapsed = prepared.preprocessing_seconds + combination_seconds
    return BoundedRepresentabilityResult(
        representables=representables,
        state_counts=state_counts,
        preprocessing_seconds=prepared.preprocessing_seconds,
        combination_seconds=combination_seconds,
        elapsed_seconds=elapsed,
        value_counts=prepared.value_counts,
        coefficient_order=prepared.coefficient_order,
        unique_preprocessed_coeffs=prepared.unique_preprocessed_coeffs,
    )


def bounded_truants_batched(
    order: RealQuadraticOrder,
    coeffs: Sequence[Pair],
    trace_bound: int,
    reorder_by_value_count: bool = True,
    cache_repeated_coeffs: bool = True,
) -> List[Pair]:
    """Return the bounded truant set by batched dynamic programming.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    coeffs : sequence[tuple[int, int]]
        Totally positive diagonal coefficients.
    trace_bound : int
        Upper trace bound.
    reorder_by_value_count : bool, default=True
        Whether to sort coefficients by increasing bounded search-set size.
    cache_repeated_coeffs : bool, default=True
        Whether to preprocess each distinct repeated coefficient only once.

    Returns
    -------
    list[tuple[int, int]]
        All totally positive omitted targets of trace at most ``trace_bound``.
    """
    targets = enumerate_totally_positive_targets(order, trace_bound)
    represented = set(
        batched_bounded_representables(
            order,
            coeffs,
            trace_bound,
            reorder_by_value_count=reorder_by_value_count,
            cache_repeated_coeffs=cache_repeated_coeffs,
        ).representables
    )
    return [alpha for alpha in targets if alpha not in represented]


def bounded_truants_naive(
    order: RealQuadraticOrder,
    coeffs: Sequence[Pair],
    trace_bound: int,
    method: str = "mitm",
    reorder_by_value_count: bool = True,
    cache_repeated_coeffs: bool = True,
) -> List[Pair]:
    """Return the bounded truant set by solving each target separately.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    coeffs : sequence[tuple[int, int]]
        Totally positive diagonal coefficients.
    trace_bound : int
        Upper trace bound.
    method : {"mitm", "dp"}, default="mitm"
        Representability solver.
    reorder_by_value_count : bool, default=True
        Whether to sort coefficients by increasing single-target search-set
        size before the combination phase.
    cache_repeated_coeffs : bool, default=True
        Whether to preprocess each distinct repeated coefficient only once.

    Returns
    -------
    list[tuple[int, int]]
        All omitted bounded targets.
    """
    targets = enumerate_totally_positive_targets(order, trace_bound)
    result: List[Pair] = []
    for alpha in targets:
        if method == "mitm":
            represented = diagonal_representability_mitm(
                order,
                coeffs,
                alpha,
                reorder_by_value_count=reorder_by_value_count,
                cache_repeated_coeffs=cache_repeated_coeffs,
            ).represented
        elif method == "dp":
            represented = diagonal_representability_dp(
                order,
                coeffs,
                alpha,
                reorder_by_value_count=reorder_by_value_count,
                cache_repeated_coeffs=cache_repeated_coeffs,
            ).represented
        else:
            raise ValueError(f"Unknown method: {method}")
        if not represented:
            result.append(alpha)
    return result


def row_candidate_count_by_trace(
    order: RealQuadraticOrder,
    coeff: Pair,
    trace_bound: int,
) -> int:
    """Return the exact number of canonical row candidates of trace at most ``trace_bound``.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    coeff : tuple[int, int]
        Totally positive coefficient.
    trace_bound : int
        Trace bound.

    Returns
    -------
    int
        Number of canonical non-zero root candidates in the trace ellipse.
    """
    return enumerate_weighted_values_by_trace(order, coeff, trace_bound, return_all_roots=False).stats.trace_candidates


def sharp_box_candidate_count_by_trace(
    order: RealQuadraticOrder,
    coeff: Pair,
    trace_bound: int,
) -> int:
    """Return the canonical sharp-box candidate count for ``q_a <= trace_bound``.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    coeff : tuple[int, int]
        Totally positive coefficient.
    trace_bound : int
        Trace bound.

    Returns
    -------
    int
        Number of canonical candidates in the sharp axis-aligned box.
    """
    form = weighted_trace_form(order, coeff)
    m_max, n_max = form.sharp_box_bounds(trace_bound)
    return m_max + n_max * (2 * m_max + 1)


def balancing_benchmark_row(
    order: RealQuadraticOrder,
    coeff: Pair,
    trace_bound: int,
) -> BalancingBenchmarkRow:
    """Return a deterministic balancing-comparison row.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    coeff : tuple[int, int]
        Totally positive coefficient.
    trace_bound : int
        Trace bound for the ellipse search.

    Returns
    -------
    BalancingBenchmarkRow
        Candidate-count comparison before and after unit-square balancing.
    """
    balanced_coeff, _ = order.balance_by_unit_square(coeff)
    return BalancingBenchmarkRow(
        D=order.D,
        coeff=coeff,
        balanced_coeff=balanced_coeff,
        trace_bound=trace_bound,
        original_row_candidates=row_candidate_count_by_trace(order, coeff, trace_bound),
        balanced_row_candidates=row_candidate_count_by_trace(order, balanced_coeff, trace_bound),
        original_sharp_box_candidates=sharp_box_candidate_count_by_trace(order, coeff, trace_bound),
        balanced_sharp_box_candidates=sharp_box_candidate_count_by_trace(order, balanced_coeff, trace_bound),
    )


def benchmark_optimisation(
    cases: Iterable[Tuple[int, Sequence[Pair], Pair]],
    repeats: int = 5,
) -> List[OptimisationBenchmarkRow]:
    """Benchmark coefficient reordering and repeated-coefficient caching.

    Parameters
    ----------
    cases : iterable[tuple[int, sequence[tuple[int, int]], tuple[int, int]]]
        Triples ``(D, coeffs, alpha)``.
    repeats : int, default=5
        Number of timing repetitions.

    Returns
    -------
    list[OptimisationBenchmarkRow]
        Optimisation benchmark rows.
    """
    rows: List[OptimisationBenchmarkRow] = []
    for D, coeffs, alpha in cases:
        order = RealQuadraticOrder(D)
        original = diagonal_representability_dp(
            order,
            coeffs,
            alpha,
            reorder_by_value_count=False,
            cache_repeated_coeffs=False,
        )
        optimised = diagonal_representability_dp(
            order,
            coeffs,
            alpha,
            reorder_by_value_count=True,
            cache_repeated_coeffs=True,
        )
        uncached_pre = [
            _prepare_single_target(
                order,
                coeffs,
                alpha,
                reorder_by_value_count=False,
                cache_repeated_coeffs=False,
            ).preprocessing_seconds
            for _ in range(repeats)
        ]
        cached_pre = [
            _prepare_single_target(
                order,
                coeffs,
                alpha,
                reorder_by_value_count=True,
                cache_repeated_coeffs=True,
            ).preprocessing_seconds
            for _ in range(repeats)
        ]
        original_total = [
            diagonal_representability_dp(
                order,
                coeffs,
                alpha,
                reorder_by_value_count=False,
                cache_repeated_coeffs=False,
            ).elapsed_seconds
            for _ in range(repeats)
        ]
        optimised_total = [
            diagonal_representability_dp(
                order,
                coeffs,
                alpha,
                reorder_by_value_count=True,
                cache_repeated_coeffs=True,
            ).elapsed_seconds
            for _ in range(repeats)
        ]
        rows.append(
            OptimisationBenchmarkRow(
                D=D,
                coeffs=tuple(coeffs),
                alpha=alpha,
                original_distinct_value_counts=original.distinct_value_counts,
                optimised_distinct_value_counts=optimised.distinct_value_counts,
                original_state_counts=original.state_counts,
                optimised_state_counts=optimised.state_counts,
                coefficient_order=optimised.coefficient_order,
                uncached_preprocessing_ms=1000.0 * median(uncached_pre),
                cached_preprocessing_ms=1000.0 * median(cached_pre),
                original_total_ms=1000.0 * median(original_total),
                optimised_total_ms=1000.0 * median(optimised_total),
            )
        )
    return rows


def benchmark_weighted_enumeration(
    cases: Iterable[Tuple[int, Pair, Pair]],
    repeats: int = 5,
) -> List[EnumerationBenchmarkRow]:
    """Benchmark weighted enumeration strategies.

    Parameters
    ----------
    cases : iterable[tuple[int, tuple[int, int], tuple[int, int]]]
        Triples ``(D, coeff, alpha)``.
    repeats : int, default=5
        Number of timing repetitions.

    Returns
    -------
    list[EnumerationBenchmarkRow]
        Enumeration benchmark rows.
    """
    rows: List[EnumerationBenchmarkRow] = []
    for D, coeff, alpha in cases:
        order = RealQuadraticOrder(D)
        row = enumerate_weighted_search(order, coeff, alpha, return_all_roots=False)
        sharp = enumerate_weighted_search_box(order, coeff, alpha, use_sharp_box=True, return_all_roots=False)
        crude = enumerate_weighted_search_box(order, coeff, alpha, use_sharp_box=False, return_all_roots=False)
        if set(row.value_to_root) != set(sharp.value_to_root) or set(row.value_to_root) != set(crude.value_to_root):
            raise AssertionError("Enumeration methods disagree on the weighted values")
        row_times = [enumerate_weighted_search(order, coeff, alpha, return_all_roots=False).stats.elapsed_seconds for _ in range(repeats)]
        sharp_times = [
            enumerate_weighted_search_box(order, coeff, alpha, use_sharp_box=True, return_all_roots=False).stats.elapsed_seconds
            for _ in range(repeats)
        ]
        crude_times = [
            enumerate_weighted_search_box(order, coeff, alpha, use_sharp_box=False, return_all_roots=False).stats.elapsed_seconds
            for _ in range(repeats)
        ]
        rows.append(
            EnumerationBenchmarkRow(
                D=D,
                coeff=coeff,
                alpha=alpha,
                distinct_values=row.stats.distinct_values,
                accepted_roots=row.stats.accepted_roots,
                sharp_box_candidates=sharp.stats.trace_candidates,
                crude_box_candidates=crude.stats.trace_candidates,
                row_candidates=row.stats.trace_candidates,
                sharp_box_ms=1000.0 * median(sharp_times),
                crude_box_ms=1000.0 * median(crude_times),
                row_ms=1000.0 * median(row_times),
            )
        )
    return rows


def benchmark_representability(
    cases: Iterable[Tuple[int, Sequence[Pair], Pair]],
    repeats: int = 5,
) -> List[RepresentabilityBenchmarkRow]:
    """Benchmark the DP and meet-in-the-middle representability layers.

    Parameters
    ----------
    cases : iterable[tuple[int, sequence[tuple[int, int]], tuple[int, int]]]
        Triples ``(D, coeffs, alpha)``.
    repeats : int, default=5
        Number of timing repetitions.

    Returns
    -------
    list[RepresentabilityBenchmarkRow]
        Benchmark rows for the representability layer.
    """
    rows: List[RepresentabilityBenchmarkRow] = []
    for D, coeffs, alpha in cases:
        order = RealQuadraticOrder(D)
        dp_result = diagonal_representability_dp(order, coeffs, alpha)
        mitm_result = diagonal_representability_mitm(order, coeffs, alpha)
        if dp_result.represented != mitm_result.represented:
            raise AssertionError("DP and MITM disagree on representability")

        prepared = _prepare_single_target(order, coeffs, alpha, reorder_by_value_count=True, cache_repeated_coeffs=True)
        represented_dp, _, dp_state_counts, _ = _dp_from_value_lists(order, prepared.value_lists, alpha)
        left_values = prepared.value_lists[: len(coeffs) // 2]
        right_values = prepared.value_lists[len(coeffs) // 2 :]
        left_sums, _ = _half_sums(order, left_values, alpha)
        right_sums, _ = _half_sums(order, right_values, alpha)
        represented_mitm = False
        if len(left_sums) <= len(right_sums):
            small_sums = left_sums
            large_sums = right_sums
        else:
            small_sums = right_sums
            large_sums = left_sums
        for partial_sum in small_sums:
            if order.sub(alpha, partial_sum) in large_sums:
                represented_mitm = True
                break
        if represented_dp != dp_result.represented or represented_mitm != mitm_result.represented:
            raise AssertionError("Combination-only routines disagree with the public API")

        dp_runs = [diagonal_representability_dp(order, coeffs, alpha) for _ in range(repeats)]
        mitm_runs = [diagonal_representability_mitm(order, coeffs, alpha) for _ in range(repeats)]

        rows.append(
            RepresentabilityBenchmarkRow(
                D=D,
                coeffs=tuple(coeffs),
                alpha=alpha,
                coefficient_order=prepared.coefficient_order,
                distinct_value_counts=prepared.distinct_value_counts,
                dp_state_counts=dp_state_counts,
                mitm_left_states=len(left_sums),
                mitm_right_states=len(right_sums),
                preprocessing_ms=1000.0 * median(run.preprocessing_seconds for run in dp_runs),
                dp_combination_ms=1000.0 * median(run.combination_seconds for run in dp_runs),
                dp_total_ms=1000.0 * median(run.elapsed_seconds for run in dp_runs),
                mitm_combination_ms=1000.0 * median(run.combination_seconds for run in mitm_runs),
                mitm_total_ms=1000.0 * median(run.elapsed_seconds for run in mitm_runs),
            )
        )
    return rows


def benchmark_bounded_truants(
    cases: Iterable[Tuple[int, Sequence[Pair], int]],
    repeats: int = 3,
) -> List[BoundedBenchmarkRow]:
    """Benchmark batched and naive bounded-truant computation.

    Parameters
    ----------
    cases : iterable[tuple[int, sequence[tuple[int, int]], int]]
        Triples ``(D, coeffs, B)``.
    repeats : int, default=3
        Number of timing repetitions.

    Returns
    -------
    list[BoundedBenchmarkRow]
        Benchmark rows for bounded truant search.
    """
    rows: List[BoundedBenchmarkRow] = []
    for D, coeffs, trace_bound in cases:
        order = RealQuadraticOrder(D)
        batched = batched_bounded_representables(order, coeffs, trace_bound)
        batched_truants = bounded_truants_batched(order, coeffs, trace_bound)
        naive_truants = bounded_truants_naive(order, coeffs, trace_bound, method="mitm")
        if batched_truants != naive_truants:
            raise AssertionError("Batched and naive bounded truant search disagree")

        prepared = _prepare_bounded_value_lists(order, coeffs, trace_bound, reorder_by_value_count=True, cache_repeated_coeffs=True)
        _, state_counts, _ = _batched_bounded_from_value_lists(order, prepared.value_lists, trace_bound)

        batched_preprocessing_times = [
            _prepare_bounded_value_lists(order, coeffs, trace_bound, reorder_by_value_count=True, cache_repeated_coeffs=True).preprocessing_seconds
            for _ in range(repeats)
        ]
        batched_combination_times = [
            _batched_bounded_from_value_lists(order, prepared.value_lists, trace_bound)[2]
            for _ in range(repeats)
        ]
        batched_total_times = []
        naive_total_times = []
        for _ in range(repeats):
            start = perf_counter()
            bounded_truants_batched(order, coeffs, trace_bound)
            batched_total_times.append(perf_counter() - start)
            start = perf_counter()
            bounded_truants_naive(order, coeffs, trace_bound, method="mitm", reorder_by_value_count=True, cache_repeated_coeffs=True)
            naive_total_times.append(perf_counter() - start)
        rows.append(
            BoundedBenchmarkRow(
                D=D,
                coeffs=tuple(coeffs),
                trace_bound=trace_bound,
                representables=len(batched.representables),
                truants=len(batched_truants),
                batched_state_counts=state_counts,
                value_counts=prepared.value_counts,
                batched_preprocessing_ms=1000.0 * median(batched_preprocessing_times),
                batched_combination_ms=1000.0 * median(batched_combination_times),
                batched_total_ms=1000.0 * median(batched_total_times),
                naive_total_ms=1000.0 * median(naive_total_times),
            )
        )
    return rows


def benchmark_generic_baseline(
    cases: Iterable[Tuple[int, Sequence[Pair], Pair]],
    repeats: int = 3,
) -> List[GenericBaselineBenchmarkRow]:
    """Benchmark a generic exact representability baseline.

    The benchmark shares the same exact preprocessing stage as the specialised
    representability routines and then replaces the DP / meet-in-the-middle
    layer by explicit Cartesian-product search over the weighted value lists.
    This yields a fully reproducible surrogate for a direct general-purpose CAS
    workflow while keeping the artifact self-contained.

    Parameters
    ----------
    cases : iterable[tuple[int, sequence[tuple[int, int]], tuple[int, int]]]
        Triples ``(D, coeffs, alpha)``.
    repeats : int, default=3
        Number of timing repetitions.

    Returns
    -------
    list[GenericBaselineBenchmarkRow]
        Benchmark rows for the generic exact baseline.
    """
    rows: List[GenericBaselineBenchmarkRow] = []
    for D, coeffs, alpha in cases:
        order = RealQuadraticOrder(D)
        dp_result = diagonal_representability_dp(order, coeffs, alpha)
        mitm_result = diagonal_representability_mitm(order, coeffs, alpha)
        prepared = _prepare_single_target(
            order,
            coeffs,
            alpha,
            reorder_by_value_count=True,
            cache_repeated_coeffs=True,
        )
        represented_generic, _, cartesian_product_size, _ = _cartesian_product_from_value_lists(
            order, prepared.value_lists, alpha
        )
        if (
            represented_generic != dp_result.represented
            or represented_generic != mitm_result.represented
        ):
            raise AssertionError("Generic exact baseline disagrees with DP/MITM")
        # Recompute preprocessing and full Cartesian-product search in each repetition.
        # Keep paired preprocessing timings with the total timings so the reported
        # generic total is not accidentally smaller than an independently sampled
        # preprocessing median on very small benchmark instances.
        preprocessing_times: List[float] = []
        generic_total_times: List[float] = []
        for _ in range(repeats):
            prepared_rt = _prepare_single_target(
                order,
                coeffs,
                alpha,
                reorder_by_value_count=True,
                cache_repeated_coeffs=True,
            )
            _, _, _, combination_rt = _cartesian_product_from_value_lists(
                order, prepared_rt.value_lists, alpha
            )
            preprocessing_times.append(prepared_rt.preprocessing_seconds)
            generic_total_times.append(prepared_rt.preprocessing_seconds + combination_rt)
        dp_total_times = [diagonal_representability_dp(order, coeffs, alpha).elapsed_seconds for _ in range(repeats)]
        mitm_total_times = [diagonal_representability_mitm(order, coeffs, alpha).elapsed_seconds for _ in range(repeats)]
        rows.append(
            GenericBaselineBenchmarkRow(
                D=D,
                coeffs=tuple(coeffs),
                alpha=alpha,
                coefficient_order=prepared.coefficient_order,
                distinct_value_counts=prepared.distinct_value_counts,
                cartesian_product_size=cartesian_product_size,
                preprocessing_ms=1000.0 * median(preprocessing_times),
                dp_total_ms=1000.0 * median(dp_total_times),
                mitm_total_ms=1000.0 * median(mitm_total_times),
                generic_total_ms=1000.0 * median(generic_total_times),
            )
        )
    return rows


def format_pair(x: Pair) -> str:
    """Return a compact string representation of a coefficient pair.

    Parameters
    ----------
    x : tuple[int, int]
        Coefficient pair.

    Returns
    -------
    str
        Human-readable string representation.
    """
    return f"({x[0]},{x[1]})"


def _is_squarefree(n: int) -> bool:
    """Return whether ``n`` is squarefree.

    Parameters
    ----------
    n : int
        Positive integer.

    Returns
    -------
    bool
        Whether no prime square divides ``n``.
    """
    if n < 1:
        return False
    if n % 4 == 0:
        return False
    m = n
    p = 2
    while p * p <= m:
        exponent = 0
        while m % p == 0:
            exponent += 1
            m //= p
            if exponent >= 2:
                return False
        p = 3 if p == 2 else p + 2
    return True


def _minimal_pell_solution(D: int) -> Tuple[int, int]:
    """Return the minimal positive solution of ``x^2 - D y^2 = 1``.

    Parameters
    ----------
    D : int
        Positive non-square integer.

    Returns
    -------
    tuple[int, int]
        The minimal positive solution.
    """
    a0 = math.isqrt(D)
    if a0 * a0 == D:
        raise ValueError("Pell equation requires non-square D")
    m = 0
    d = 1
    a = a0
    p_nm2, p_nm1 = 0, 1
    q_nm2, q_nm1 = 1, 0
    while True:
        p_n = a * p_nm1 + p_nm2
        q_n = a * q_nm1 + q_nm2
        if p_n * p_n - D * q_n * q_n == 1:
            return (p_n, q_n)
        m = d * a - m
        d = (D - m * m) // d
        a = (a0 + m) // d
        p_nm2, p_nm1 = p_nm1, p_n
        q_nm2, q_nm1 = q_nm1, q_n



def _sign_of_p_plus_q_sqrtD(p: int, q: int, D: int) -> int:
    """Return the sign of ``p + q * sqrt(D)`` for squarefree ``D >= 2``.

    Parameters
    ----------
    p, q : int
        Integer coefficients.
    D : int
        Positive non-square radicand.

    Returns
    -------
    int
        ``-1`` for negative, ``0`` for zero, and ``1`` for positive.
    """
    if p == 0 and q == 0:
        return 0
    if q == 0:
        return 1 if p > 0 else -1
    if p == 0:
        return 1 if q > 0 else -1
    if p > 0 and q > 0:
        return 1
    if p < 0 and q < 0:
        return -1
    lhs = p * p
    rhs = D * q * q
    if lhs == rhs:
        return 0
    if q > 0:
        # Here p < 0.
        return 1 if rhs > lhs else -1
    # Here q < 0 and therefore p > 0.
    return 1 if lhs > rhs else -1


def _floor_div(n: int, d: int) -> int:
    """Return ``floor(n / d)`` for ``d > 0``.

    Parameters
    ----------
    n : int
        Numerator.
    d : int
        Positive denominator.

    Returns
    -------
    int
        The floor quotient.
    """
    return n // d


def _ceil_div(n: int, d: int) -> int:
    """Return ``ceil(n / d)`` for ``d > 0``.

    Parameters
    ----------
    n : int
        Numerator.
    d : int
        Positive denominator.

    Returns
    -------
    int
        The ceiling quotient.
    """
    return -((-n) // d)
