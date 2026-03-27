"""Pythagoras number P(O_D) and explicit Pythagoras elements (Algorithm 4).

The classification is from [KRS, Theorem 3.1].  Witnesses are verified
via the exact_length routine.
"""

from __future__ import annotations

from typing import Tuple

from quadratic_sos.ring import RealQuadraticOrder
from quadratic_sos.length import exact_length

Pair = Tuple[int, int]


def pythagoras_number(order: RealQuadraticOrder) -> int:
    """Return P(O_D)."""
    return order.pythagoras_number()


def pythagoras_element(order: RealQuadraticOrder) -> Pair:
    """Return a coefficient pair (a, b) for an element a_D with ℓ_D(a_D) = P(O_D)."""
    D = order.D

    if D == 2:
        return (6, 2)
    if D == 3:
        return (9, 4)
    if D == 5:
        # (7+√5)/2 = 3 + 1·ω₅ in the basis {1, ω₅}
        return (3, 1)
    if D == 6:
        return (10, 2)
    if D == 7:
        return (11, 2)

    if D == 13:
        # 12 + 2√13.  In the basis {1, ω₁₃=(1+√13)/2}:
        # 12 + 2√13 = 12 + 2·(2ω₁₃ - 1) = 10 + 4ω₁₃
        return (10, 4)

    # Generic P = 5 witness
    if order.is_mod1:
        # 7 + ω_D² = 7 + (D-1)/4 + ω_D = (7 + (D-1)//4, 1)
        # since ω² = (D-1)/4 + ω
        return (7 + (D - 1) // 4, 1)
    else:
        # 7 + (1+√D)² = 7 + 1 + D + 2√D = (8+D, 2)
        return (8 + D, 2)


def verify_pythagoras_element(order: RealQuadraticOrder) -> Tuple[int, Pair, int]:
    """Return (P, witness, computed_length) as a self-check.

    The computed_length should equal P.  If it doesn't, something is wrong.
    """
    p = pythagoras_number(order)
    w = pythagoras_element(order)
    ell = exact_length(order, w)
    return (p, w, ell)
