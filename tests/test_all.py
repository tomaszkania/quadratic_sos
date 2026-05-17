"""Regression tests for the ``quadratic_sos`` package.

The tests reproduce the paper's worked examples, witness elements, validation
sweep, and the ``D=10`` distribution table.
"""

from __future__ import annotations

from typing import Callable

import pytest

from quadratic_sos import RealQuadraticOrder
from quadratic_sos.enumeration import enumerate_search_set, is_representable, relevant_squares
from quadratic_sos.experiments import search_profile, validation_sweep
from quadratic_sos.length import decomposition, exact_length
from quadratic_sos.pythagoras import pythagoras_number, verify_pythagoras_element
from quadratic_sos.tables import length_distribution


def test_norm_trace_D10() -> None:
    """Check norm and trace in ``O_10``."""
    order = RealQuadraticOrder(10)
    assert order.norm((18, 2)) == 284
    assert order.trace((18, 2)) == 36


def test_norm_trace_D5() -> None:
    """Check norm and trace in the integral basis of ``O_5``."""
    order = RealQuadraticOrder(5)
    assert order.norm((3, 1)) == 11
    assert order.trace((3, 1)) == 7


def test_square_D10() -> None:
    """Check squaring in ``Z[sqrt(10)]``."""
    order = RealQuadraticOrder(10)
    assert order.sqr((1, 1)) == (11, 2)


def test_square_D5() -> None:
    """Check squaring in ``Z[(1+sqrt(5))/2]``."""
    order = RealQuadraticOrder(5)
    assert order.sqr((0, 1)) == (1, 1)


def test_squarefree_validation() -> None:
    """Reject non-squarefree radicands."""
    with pytest.raises(ValueError, match="not squarefree"):
        RealQuadraticOrder(12)


def test_mixed_order_arithmetic_rejected() -> None:
    """Reject arithmetic between elements from different orders."""
    x = RealQuadraticOrder(5).elem(1, 0)
    y = RealQuadraticOrder(10).elem(1, 0)
    with pytest.raises(ValueError, match="same real quadratic order"):
        _ = x + y


def test_immediate_obstruction_D10() -> None:
    """Example 7.1: ``3 + sqrt(10)`` is not representable."""
    order = RealQuadraticOrder(10)
    assert not is_representable(order, (3, 1))


def test_representable_D10() -> None:
    """Check that the length-five example is representable."""
    order = RealQuadraticOrder(10)
    assert is_representable(order, (18, 2))


def test_not_totally_positive() -> None:
    """Non-totally-positive elements are not sums of squares."""
    order = RealQuadraticOrder(10)
    assert not is_representable(order, (-1, 0))


def test_search_set_D10() -> None:
    """Reproduce ``S(18 + 2 sqrt(10))`` from the worked example."""
    order = RealQuadraticOrder(10)
    alpha = (18, 2)
    roots = enumerate_search_set(order, alpha)
    squares = relevant_squares(order, alpha)
    expected_squares = {(0, 0), (1, 0), (4, 0), (9, 0), (10, 0), (11, 2)}
    assert len(roots) == 11
    assert set(squares) == expected_squares


def test_length_1() -> None:
    """The element ``1`` has exact length one."""
    order = RealQuadraticOrder(10)
    assert exact_length(order, (1, 0)) == 1


def test_length_2() -> None:
    """The element ``2`` has exact length two."""
    order = RealQuadraticOrder(10)
    assert exact_length(order, (2, 0)) == 2


def test_length_3() -> None:
    """The element ``3`` has exact length three."""
    order = RealQuadraticOrder(10)
    assert exact_length(order, (3, 0)) == 3


def test_length_4() -> None:
    """The element ``7`` has exact length four."""
    order = RealQuadraticOrder(10)
    assert exact_length(order, (7, 0)) == 4


def test_length_5_D10() -> None:
    """The element ``18 + 2 sqrt(10)`` has exact length five."""
    order = RealQuadraticOrder(10)
    assert exact_length(order, (18, 2)) == 5


def test_length_infinity() -> None:
    """Non-representability is returned as ``None``."""
    order = RealQuadraticOrder(10)
    assert exact_length(order, (3, 1)) is None


def test_length_zero_rejected() -> None:
    """The zero element is outside the package's exact-length API."""
    order = RealQuadraticOrder(10)
    with pytest.raises(ValueError, match="nonzero"):
        exact_length(order, (0, 0))
    with pytest.raises(ValueError, match="nonzero"):
        decomposition(order, (0, 0))


def test_length_D2() -> None:
    """Check a ternary Pythagoras witness."""
    order = RealQuadraticOrder(2)
    assert exact_length(order, (6, 2)) == 3


def test_length_D6() -> None:
    """Check a quaternary Pythagoras witness."""
    order = RealQuadraticOrder(6)
    assert exact_length(order, (10, 2)) == 4


def test_decomposition_verifies() -> None:
    """Recovered decompositions sum back to their target element."""
    order = RealQuadraticOrder(10)
    for alpha in [(1, 0), (2, 0), (3, 0), (7, 0)]:
        recovered = decomposition(order, alpha)
        assert recovered is not None
        total = (0, 0)
        for root in recovered:
            total = order.add(total, order.sqr(root))
        assert total == alpha


def test_pythagoras_numbers() -> None:
    """Reproduce the real-quadratic Pythagoras-number classification."""
    cases = {2: 3, 3: 3, 5: 3, 6: 4, 7: 4, 10: 5, 11: 5, 13: 5, 17: 5}
    for radicand, expected in cases.items():
        order = RealQuadraticOrder(radicand)
        assert pythagoras_number(order) == expected


def test_pythagoras_witnesses() -> None:
    """Verify that every bundled witness has maximal exact length."""
    for radicand in [2, 3, 5, 6, 7, 10, 11, 13, 17, 21, 29]:
        order = RealQuadraticOrder(radicand)
        pythagoras, _, length = verify_pythagoras_element(order)
        assert length == pythagoras


def test_unbalanced_length_matches_balanced() -> None:
    """Unit-square balancing preserves exact length."""
    order = RealQuadraticOrder(10)
    beta = (18, 2)
    unit_square = order.sqr(order.pell_unit)
    alpha = order.mul(order.mul(beta, unit_square), unit_square)
    assert exact_length(order, alpha) == exact_length(order, alpha, balance=False) == 5


def test_balancing_reduces_trace_row_bound() -> None:
    """Balancing shrinks the trace-row bound on a Pell orbit."""
    order = RealQuadraticOrder(10)
    beta = (18, 2)
    unit_square = order.sqr(order.pell_unit)
    alpha = order.mul(order.mul(order.mul(beta, unit_square), unit_square), unit_square)
    raw = search_profile(order, alpha, balance=False)
    balanced = search_profile(order, alpha, balance=True)
    assert raw.trace_v_bound > balanced.trace_v_bound
    assert raw.box_size == balanced.box_size == 11


def test_validation_sweep_small() -> None:
    """Small brute-force validation sweep has no mismatches."""
    summary = validation_sweep(13, 12)
    assert summary.mismatches == 0
    assert summary.field_count == 8
    assert summary.element_count > 0


def test_table_D10_trace40() -> None:
    """Reproduce the ``D=10`` trace-40 distribution table."""
    order = RealQuadraticOrder(10)
    counts, _ = length_distribution(order, 40)
    assert sum(counts.values()) == 134
    assert counts[1] == 11
    assert counts[2] == 22
    assert counts[3] == 14
    assert counts[4] == 9
    assert counts[5] == 2
    assert counts[None] == 76


if __name__ == "__main__":
    test_functions: list[Callable[[], None]] = [
        value for key, value in sorted(globals().items()) if key.startswith("test_")
    ]
    passed = 0
    for function in test_functions:
        try:
            function()
            print(f"  ✓ {function.__name__}")
            passed += 1
        except Exception as exc:  # pragma: no cover - manual test runner
            print(f"  ✗ {function.__name__}: {exc}")
    print(f"\n{passed}/{len(test_functions)} tests passed.")
