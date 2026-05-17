"""Exact arithmetic in the ring of integers ``O_D`` of ``Q(sqrt(D))``.

All computations use Python's arbitrary-precision integers. Elements are
represented by coefficient pairs ``(a, b)`` relative to the integral basis
``{1, omega_D}``. The public API is restricted to squarefree radicands ``D``
and to maximal real quadratic orders.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import cached_property
from typing import Tuple

Pair = Tuple[int, int]


@dataclass(frozen=True)
class RealQuadraticOrder:
    """The maximal order ``O_D = Z[omega_D]`` of ``Q(sqrt(D))``.

    Parameters
    ----------
    D : int
        Squarefree integer at least ``2``.
    """

    D: int

    def __post_init__(self) -> None:
        if self.D < 2:
            raise ValueError(f"D must be >= 2, got {self.D}")
        s = math.isqrt(self.D)
        if s * s == self.D:
            raise ValueError(f"D = {self.D} is a perfect square")
        if not _is_squarefree(self.D):
            raise ValueError(f"D = {self.D} is not squarefree")

    @property
    def is_mod1(self) -> bool:
        """Return ``True`` when ``D ≡ 1 (mod 4)``."""
        return self.D % 4 == 1

    @cached_property
    def omega_minpoly(self) -> Tuple[int, int]:
        """Return ``(p, q)`` such that ``omega_D^2 = p + q * omega_D``."""
        if self.is_mod1:
            return ((self.D - 1) // 4, 1)
        return (self.D, 0)

    @property
    def one(self) -> Pair:
        """Return the multiplicative identity as a coefficient pair."""
        return (1, 0)

    def elem(self, a: int, b: int) -> "QuadraticInteger":
        """Create the element ``a + b * omega_D`` in ``O_D``."""
        return QuadraticInteger(a, b, self)

    def from_rational(self, n: int) -> "QuadraticInteger":
        """Create the rational integer ``n`` inside ``O_D``."""
        return QuadraticInteger(n, 0, self)

    def add(self, x: Pair, y: Pair) -> Pair:
        return (x[0] + y[0], x[1] + y[1])

    def sub(self, x: Pair, y: Pair) -> Pair:
        return (x[0] - y[0], x[1] - y[1])

    def mul(self, x: Pair, y: Pair) -> Pair:
        """Multiply two elements represented as coefficient pairs."""
        a1, b1 = x
        a2, b2 = y
        p, q = self.omega_minpoly
        return (a1 * a2 + b1 * b2 * p, a1 * b2 + a2 * b1 + b1 * b2 * q)

    def sqr(self, x: Pair) -> Pair:
        """Square an element represented as a coefficient pair."""
        u, v = x
        p, q = self.omega_minpoly
        return (u * u + v * v * p, 2 * u * v + v * v * q)

    def pow(self, x: Pair, n: int) -> Pair:
        """Return ``x**n`` for ``n >= 0`` using binary exponentiation.

        Parameters
        ----------
        x : tuple[int, int]
            Base element in coefficient form.
        n : int
            Non-negative exponent.

        Returns
        -------
        tuple[int, int]
            The power ``x**n``.
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

    def norm(self, x: Pair) -> int:
        """Return the field norm ``N(x)`` as an integer."""
        a, b = x
        if self.is_mod1:
            return a * a + a * b - ((self.D - 1) // 4) * b * b
        return a * a - self.D * b * b

    def trace(self, x: Pair) -> int:
        """Return the field trace ``Tr(x)`` as an integer."""
        a, b = x
        return 2 * a + b if self.is_mod1 else 2 * a

    def conjugate(self, x: Pair) -> Pair:
        """Return the Galois conjugate of ``x`` in coefficient form."""
        a, b = x
        if self.is_mod1:
            return (a + b, -b)
        return (a, -b)

    def _sigma1_sign(self, x: Pair) -> int:
        """Return the sign of the first real embedding ``sigma_1(x)``.

        Returns
        -------
        int
            ``-1``, ``0`` or ``1``.
        """
        a, b = x
        if self.is_mod1:
            r, s = 2 * a + b, b
        else:
            r, s = a, b
        return _sign_of_r_plus_s_sqrtD(r, s, self.D)

    def _sigma2_sign(self, x: Pair) -> int:
        """Return the sign of the second real embedding ``sigma_2(x)``."""
        return self._sigma1_sign(self.conjugate(x))

    def is_totally_positive(self, x: Pair) -> bool:
        """Return ``True`` when both real embeddings of ``x`` are positive."""
        return self._sigma1_sign(x) > 0 and self._sigma2_sign(x) > 0

    def sigma1_le(self, x: Pair, y: Pair) -> bool:
        """Return ``True`` when ``sigma_1(x) <= sigma_1(y)``."""
        return self._sigma1_sign(self.sub(y, x)) >= 0

    def sigma2_le(self, x: Pair, y: Pair) -> bool:
        """Return ``True`` when ``sigma_2(x) <= sigma_2(y)``."""
        return self._sigma2_sign(self.sub(y, x)) >= 0

    def trace_row_bounds(self, alpha: Pair, v: int) -> Tuple[int, int] | None:
        """Return an exact trace-based enclosure for admissible ``u`` on row ``v``.

        If ``x = u + v*omega_D`` satisfies ``x^2 <= alpha`` in both embeddings,
        then ``Tr(x^2) <= Tr(alpha)``. This yields a finite interval of possible
        ``u`` values that can be computed exactly.

        Parameters
        ----------
        alpha : tuple[int, int]
            Totally positive target element.
        v : int
            Fixed row index.

        Returns
        -------
        tuple[int, int] or None
            A guaranteed enclosing interval ``[u_min, u_max]`` for admissible
            integers ``u`` on the row ``v``. Returns ``None`` when the trace
            inequality already forbids all solutions.
        """
        tr_alpha = self.trace(alpha)
        if self.is_mod1:
            rhs = 2 * tr_alpha - self.D * v * v
            if rhs < 0:
                return None
            s = math.isqrt(rhs)
            return (_ceil_div(-v - s, 2), _floor_div(-v + s, 2))

        rhs = tr_alpha // 2 - self.D * v * v
        if rhs < 0:
            return None
        s = math.isqrt(rhs)
        return (-s, s)

    def trace_v_bound(self, alpha: Pair) -> int:
        """Return an exact trace-based bound for ``|v|`` in the search set."""
        tr_alpha = self.trace(alpha)
        if self.is_mod1:
            return math.isqrt((2 * tr_alpha) // self.D)
        return math.isqrt((tr_alpha // 2) // self.D)

    @cached_property
    def pell_unit(self) -> Pair:
        """Return a positive norm-one unit obtained from Pell's equation.

        The returned unit is the element corresponding to the minimal positive
        solution of ``x^2 - D y^2 = 1``. For ``D ≡ 1 (mod 4)`` this need not be
        the fundamental unit of ``O_D``, but it is still a positive norm-one
        unit and therefore suffices for exact unit-square balancing.

        Notes
        -----
        The continued-fraction computation of the Pell solution is fast for the
        moderate discriminants used in the paper and test suite. For extremely
        large ``D`` it can become a noticeable setup cost before the actual
        search begins.
        """
        x, y = _minimal_pell_solution(self.D)
        if self.is_mod1:
            return (x - y, 2 * y)
        return (x, y)

    def balance_by_pell_unit(self, alpha: Pair) -> Tuple[Pair, Pair]:
        """Balance a totally positive element by even powers of a norm-one unit.

        Parameters
        ----------
        alpha : tuple[int, int]
            Totally positive input element.

        Returns
        -------
        tuple[tuple[int, int], tuple[int, int]]
            A pair ``(beta, scale)`` such that ``beta = alpha * u^(2k)`` for a
            suitable integer ``k`` and the Pell unit ``u = self.pell_unit``,
            while ``scale = u^(-k)``. Consequently, every decomposition
            ``beta = y_1^2 + ... + y_n^2`` yields
            ``alpha = (scale*y_1)^2 + ... + (scale*y_n)^2``.
        """
        if not self.is_totally_positive(alpha):
            return alpha, self.one

        unit = self.pell_unit
        unit_inv = self.conjugate(unit)
        unit_sq = self.sqr(unit)
        unit_sq_inv = self.sqr(unit_inv)

        beta = alpha
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

    def _ratio_exceeds_square(self, beta: Pair, unit_sq: Pair) -> bool:
        """Return ``True`` when ``sigma_1(beta) / sigma_2(beta)`` is too large."""
        rhs = self.mul(unit_sq, self.conjugate(beta))
        return not self.sigma1_le(beta, rhs)

    def _inverse_ratio_exceeds_square(self, beta: Pair, unit_sq: Pair) -> bool:
        """Return ``True`` when ``sigma_2(beta) / sigma_1(beta)`` is too large."""
        rhs = self.mul(unit_sq, beta)
        return not self.sigma1_le(self.conjugate(beta), rhs)

    def pythagoras_number(self) -> int:
        """Return the integral Pythagoras number ``P(O_D)``."""
        if self.D in (2, 3, 5):
            return 3
        if self.D in (6, 7):
            return 4
        return 5


@dataclass(frozen=True)
class QuadraticInteger:
    """An element ``a + b*omega_D`` of the ring ``O_D``.

    Attributes
    ----------
    a, b : int
        Coefficients in the integral basis ``{1, omega_D}``.
    order : RealQuadraticOrder
        Ambient maximal order.
    """

    a: int
    b: int
    order: RealQuadraticOrder

    @property
    def pair(self) -> Pair:
        """Return the coefficient pair of the element."""
        return (self.a, self.b)

    def _coerce_other(self, other: object) -> "QuadraticInteger":
        """Return ``other`` after validating that the ambient orders agree."""
        if not isinstance(other, QuadraticInteger):
            raise TypeError("Arithmetic is only defined between QuadraticInteger instances")
        if self.order.D != other.order.D:
            raise ValueError("QuadraticInteger operands must belong to the same real quadratic order")
        return other

    def __add__(self, other: "QuadraticInteger") -> "QuadraticInteger":
        other_q = self._coerce_other(other)
        p = self.order.add(self.pair, other_q.pair)
        return QuadraticInteger(p[0], p[1], self.order)

    def __sub__(self, other: "QuadraticInteger") -> "QuadraticInteger":
        other_q = self._coerce_other(other)
        p = self.order.sub(self.pair, other_q.pair)
        return QuadraticInteger(p[0], p[1], self.order)

    def __mul__(self, other: "QuadraticInteger") -> "QuadraticInteger":
        other_q = self._coerce_other(other)
        p = self.order.mul(self.pair, other_q.pair)
        return QuadraticInteger(p[0], p[1], self.order)

    def __neg__(self) -> "QuadraticInteger":
        return QuadraticInteger(-self.a, -self.b, self.order)

    def square(self) -> "QuadraticInteger":
        p = self.order.sqr(self.pair)
        return QuadraticInteger(p[0], p[1], self.order)

    def norm(self) -> int:
        return self.order.norm(self.pair)

    def trace(self) -> int:
        return self.order.trace(self.pair)

    def conjugate(self) -> "QuadraticInteger":
        p = self.order.conjugate(self.pair)
        return QuadraticInteger(p[0], p[1], self.order)

    def is_totally_positive(self) -> bool:
        return self.order.is_totally_positive(self.pair)

    def __repr__(self) -> str:
        D = self.order.D
        basis = f"ω_{D}" if self.order.is_mod1 else f"√{D}"
        if self.b == 0:
            return str(self.a)
        if self.a == 0:
            if self.b == 1:
                return basis
            if self.b == -1:
                return f"-{basis}"
            return f"{self.b}·{basis}"
        sign = " + " if self.b > 0 else " - "
        coeff = abs(self.b)
        c = "" if coeff == 1 else f"{coeff}·"
        return f"{self.a}{sign}{c}{basis}"

    def __eq__(self, other: object) -> bool:
        if isinstance(other, QuadraticInteger):
            return self.a == other.a and self.b == other.b and self.order.D == other.order.D
        return NotImplemented

    def __hash__(self) -> int:
        return hash((self.a, self.b, self.order.D))


def _sign_of_r_plus_s_sqrtD(r: int, s: int, D: int) -> int:
    """Return the sign of ``r + s*sqrt(D)`` exactly.

    Parameters
    ----------
    r, s : int
        Rational coefficients.
    D : int
        Squarefree radicand.

    Returns
    -------
    int
        ``-1``, ``0`` or ``1`` according to the sign of ``r + s*sqrt(D)``.
    """
    if s == 0:
        return (r > 0) - (r < 0)
    if r == 0:
        return (s > 0) - (s < 0)
    if r > 0 and s > 0:
        return 1
    if r < 0 and s < 0:
        return -1

    r2 = r * r
    s2D = s * s * D
    if r > 0:
        if r2 > s2D:
            return 1
        if r2 < s2D:
            return -1
        return 0
    if s2D > r2:
        return 1
    if s2D < r2:
        return -1
    return 0


def _is_squarefree(n: int) -> bool:
    """Return ``True`` when ``n`` is squarefree.

    Parameters
    ----------
    n : int
        Positive integer to test.

    Returns
    -------
    bool
        ``True`` if and only if no prime square divides ``n``.
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
            m //= p
            exponent += 1
            if exponent >= 2:
                return False
        p = 3 if p == 2 else p + 2
    return True


def _minimal_pell_solution(D: int) -> Tuple[int, int]:
    """Return the minimal positive solution of ``x^2 - D y^2 = 1``.

    Parameters
    ----------
    D : int
        Non-square positive integer.

    Returns
    -------
    tuple[int, int]
        The minimal solution ``(x, y)`` with ``x > 0`` and ``y > 0``.
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


def _floor_div(n: int, d: int) -> int:
    """Return ``floor(n / d)`` for integers with ``d > 0``."""
    return n // d


def _ceil_div(n: int, d: int) -> int:
    """Return ``ceil(n / d)`` for integers with ``d > 0``."""
    return -((-n) // d)
