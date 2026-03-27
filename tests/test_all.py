"""Tests for quadratic_sos — reproducing the paper's examples and table."""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from quadratic_sos.ring import RealQuadraticOrder, QuadraticInteger
from quadratic_sos.enumeration import is_representable, relevant_squares, enumerate_search_set
from quadratic_sos.length import exact_length, decomposition
from quadratic_sos.pythagoras import pythagoras_number, pythagoras_element, verify_pythagoras_element
from quadratic_sos.tables import enumerate_totally_positive, length_distribution
from quadratic_sos.experiments import search_profile, validation_sweep


# ── Ring arithmetic ───────────────────────────────────────────────

def test_norm_trace_D10():
    O = RealQuadraticOrder(10)
    # 18 + 2√10: norm = 18² - 10·4 = 324 - 40 = 284
    assert O.norm((18, 2)) == 284
    assert O.trace((18, 2)) == 36

def test_norm_trace_D5():
    O = RealQuadraticOrder(5)
    # (7+√5)/2 = 3 + ω₅: norm = 9 + 3 - 1 = 11, trace = 7
    assert O.norm((3, 1)) == 3*3 + 3*1 - 1*1  # = 11
    assert O.trace((3, 1)) == 7

def test_square_D10():
    O = RealQuadraticOrder(10)
    # (1+√10)² = 1 + 10 + 2√10 = (11, 2)
    assert O.sqr((1, 1)) == (11, 2)

def test_square_D5():
    O = RealQuadraticOrder(5)
    # ω₅² = (D-1)/4 + ω₅ = 1 + ω₅ = (1, 1)
    assert O.sqr((0, 1)) == (1, 1)


def test_squarefree_validation():
    try:
        RealQuadraticOrder(12)
    except ValueError:
        pass
    else:
        raise AssertionError("D=12 should be rejected as non-squarefree")


def test_mixed_order_arithmetic_rejected():
    x = RealQuadraticOrder(5).elem(1, 0)
    y = RealQuadraticOrder(10).elem(1, 0)
    try:
        _ = x + y
    except ValueError:
        pass
    else:
        raise AssertionError("Mixed-order arithmetic should raise ValueError")


# ── Representability (Algorithm 1) ────────────────────────────────

def test_immediate_obstruction_D10():
    """Example 7.1: 3+√10 has odd √10-coefficient → not representable."""
    O = RealQuadraticOrder(10)
    assert not is_representable(O, (3, 1))

def test_representable_D10():
    O = RealQuadraticOrder(10)
    assert is_representable(O, (18, 2))

def test_not_tp():
    O = RealQuadraticOrder(10)
    assert not is_representable(O, (-1, 0))


# ── Enumeration ───────────────────────────────────────────────────

def test_search_set_D10():
    """End-to-end example: B(18+2√10) and S(18+2√10)."""
    O = RealQuadraticOrder(10)
    alpha = (18, 2)
    roots = enumerate_search_set(O, alpha)
    S = relevant_squares(O, alpha)
    # Paper claims S = {0, 1, 4, 9, 10, 11+2√10}
    expected_S = {(0,0), (1,0), (4,0), (9,0), (10,0), (11,2)}
    assert set(S) == expected_S


# ── Exact length (Algorithm 3) ────────────────────────────────────

def test_length_1():
    O = RealQuadraticOrder(10)
    assert exact_length(O, (1, 0)) == 1   # 1 = 1²

def test_length_2():
    O = RealQuadraticOrder(10)
    assert exact_length(O, (2, 0)) == 2   # 2 = 1² + 1²

def test_length_3():
    O = RealQuadraticOrder(10)
    assert exact_length(O, (3, 0)) == 3   # 3 = 1² + 1² + 1²

def test_length_4():
    O = RealQuadraticOrder(10)
    assert exact_length(O, (7, 0)) == 4   # 7 = 4 + 1 + 1 + 1

def test_length_5_D10():
    """End-to-end example: ℓ₁₀(18+2√10) = 5."""
    O = RealQuadraticOrder(10)
    assert exact_length(O, (18, 2)) == 5

def test_length_inf():
    O = RealQuadraticOrder(10)
    # 3+√10 is not representable
    assert exact_length(O, (3, 1)) is None

def test_length_zero_rejected():
    O = RealQuadraticOrder(10)
    try:
        exact_length(O, (0, 0))
    except ValueError:
        pass
    else:
        raise AssertionError("exact_length should reject the zero element")

    try:
        decomposition(O, (0, 0))
    except ValueError:
        pass
    else:
        raise AssertionError("decomposition should reject the zero element")

def test_length_D2():
    O = RealQuadraticOrder(2)
    # P(O₂) = 3, so max length is 3
    assert exact_length(O, (6, 2)) == 3

def test_length_D6():
    O = RealQuadraticOrder(6)
    assert exact_length(O, (10, 2)) == 4


# ── Decomposition ─────────────────────────────────────────────────

def test_decomposition_verifies():
    """Check that recovered decompositions actually sum to alpha."""
    O = RealQuadraticOrder(10)
    for alpha in [(1,0), (2,0), (3,0), (7,0)]:
        dec = decomposition(O, alpha)
        assert dec is not None
        total = (0, 0)
        for root in dec:
            total = O.add(total, O.sqr(root))
        assert total == alpha, f"Decomposition of {alpha} doesn't verify: {dec}"


# ── Pythagoras numbers ────────────────────────────────────────────

def test_pythagoras_numbers():
    cases = {2: 3, 3: 3, 5: 3, 6: 4, 7: 4, 10: 5, 11: 5, 13: 5, 17: 5}
    for D, expected in cases.items():
        O = RealQuadraticOrder(D)
        assert pythagoras_number(O) == expected

def test_pythagoras_witnesses():
    """Verify that every Pythagoras witness actually achieves the stated length."""
    for D in [2, 3, 5, 6, 7, 10, 11, 13, 17, 21, 29]:
        O = RealQuadraticOrder(D)
        p, w, ell = verify_pythagoras_element(O)
        assert ell == p, f"D={D}: P={p} but witness has length {ell}"




def test_unbalanced_length_matches_balanced():
    O = RealQuadraticOrder(10)
    beta = (18, 2)
    unit_sq = O.sqr(O.pell_unit)
    alpha = O.mul(O.mul(beta, unit_sq), unit_sq)
    assert exact_length(O, alpha) == exact_length(O, alpha, balance=False) == 5


def test_balancing_reduces_trace_row_bound():
    O = RealQuadraticOrder(10)
    beta = (18, 2)
    unit_sq = O.sqr(O.pell_unit)
    alpha = O.mul(O.mul(O.mul(beta, unit_sq), unit_sq), unit_sq)
    raw = search_profile(O, alpha, balance=False)
    bal = search_profile(O, alpha, balance=True)
    assert raw.trace_v_bound > bal.trace_v_bound
    assert raw.box_size == bal.box_size == 11


def test_validation_sweep_small():
    summary = validation_sweep(13, 12)
    assert summary.mismatches == 0
    assert summary.field_count == 8
    assert summary.element_count > 0


# ── Table reproduction ────────────────────────────────────────────

def test_table_D10_trace40():
    """Reproduce Table 1 from the paper."""
    O = RealQuadraticOrder(10)
    counts, first = length_distribution(O, 40)
    assert sum(counts.values()) == 134
    assert counts[1] == 11
    assert counts[2] == 22
    assert counts[3] == 14
    assert counts[4] == 9
    assert counts[5] == 2
    assert counts[None] == 76


# ── Run all ───────────────────────────────────────────────────────

if __name__ == "__main__":
    test_fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = 0
    for fn in test_fns:
        try:
            fn()
            print(f"  ✓ {fn.__name__}")
            passed += 1
        except Exception as e:
            print(f"  ✗ {fn.__name__}: {e}")
    print(f"\n{passed}/{len(test_fns)} tests passed.")
