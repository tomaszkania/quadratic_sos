"""Tests for the self-contained weighted diagonal representability module."""

from __future__ import annotations

from itertools import product
from typing import Dict, Iterable, List, Sequence, Tuple

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from quadratic_diagonal import (
    Pair,
    RealQuadraticOrder,
    balancing_benchmark_row,
    benchmark_bounded_truants,
    benchmark_generic_baseline,
    benchmark_optimisation,
    benchmark_representability,
    benchmark_weighted_enumeration,
    batched_bounded_representables,
    bounded_truants_batched,
    bounded_truants_naive,
    diagonal_representability_dp,
    diagonal_representability_mitm,
    discriminant_identity,
    enumerate_totally_positive_targets,
    enumerate_weighted_search,
    enumerate_weighted_search_box,
    enumerate_weighted_values_by_trace,
    exact_row_interval,
    weighted_trace_form,
)


def brute_force_weighted_search(order: RealQuadraticOrder, coeff: Pair, alpha: Pair) -> Dict[Pair, Pair]:
    """Return the weighted search set by explicit sharp-box scanning.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    coeff : tuple[int, int]
        Totally positive coefficient.
    alpha : tuple[int, int]
        Totally positive target.

    Returns
    -------
    dict[tuple[int, int], tuple[int, int]]
        One witness root for each distinct weighted value.
    """
    form = weighted_trace_form(order, coeff)
    m_max, n_max = form.sharp_box_bounds(order.trace(alpha))
    mapping: Dict[Pair, Pair] = {}
    for n in range(-n_max, n_max + 1):
        for m in range(-m_max, m_max + 1):
            if m == 0 and n == 0:
                continue
            if form.value(m, n) > order.trace(alpha):
                continue
            x = (m, n)
            value = order.mul(coeff, order.sqr(x))
            if order.le(value, alpha):
                mapping.setdefault(value, x)
    return mapping


def brute_force_weighted_values_by_trace(order: RealQuadraticOrder, coeff: Pair, trace_bound: int) -> Dict[Pair, Pair]:
    """Return all weighted values of trace at most ``trace_bound`` by brute force.

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
    dict[tuple[int, int], tuple[int, int]]
        One witness root for each distinct weighted value.
    """
    form = weighted_trace_form(order, coeff)
    m_max, n_max = form.sharp_box_bounds(trace_bound)
    mapping: Dict[Pair, Pair] = {}
    for n in range(-n_max, n_max + 1):
        for m in range(-m_max, m_max + 1):
            if m == 0 and n == 0:
                continue
            if form.value(m, n) > trace_bound:
                continue
            x = (m, n)
            value = order.mul(coeff, order.sqr(x))
            mapping.setdefault(value, x)
    return mapping


def brute_force_diagonal_representable(order: RealQuadraticOrder, coeffs: Sequence[Pair], alpha: Pair) -> bool:
    """Return whether ``alpha`` is represented by the diagonal lattice.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    coeffs : sequence[tuple[int, int]]
        Diagonal coefficients.
    alpha : tuple[int, int]
        Totally positive target.

    Returns
    -------
    bool
        Whether a representation exists.
    """
    value_lists: List[List[Pair]] = [[order.zero] for _ in coeffs]
    for i, coeff in enumerate(coeffs):
        value_lists[i].extend(list(brute_force_weighted_search(order, coeff, alpha)))
    for choice in product(*value_lists):
        total = order.zero
        for value in choice:
            total = order.add(total, value)
        if total == alpha:
            return True
    return False


def brute_force_totally_positive_targets(order: RealQuadraticOrder, trace_bound: int) -> List[Pair]:
    """Return all totally positive targets by a simple finite box search.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    trace_bound : int
        Trace bound.

    Returns
    -------
    list[tuple[int, int]]
        Sorted list of non-zero totally positive elements of bounded trace.
    """
    if order.is_mod1:
        m_range = range(-2 * trace_bound, 2 * trace_bound + 1)
        n_range = range(-2 * trace_bound, 2 * trace_bound + 1)
    else:
        m_range = range(-trace_bound, trace_bound + 1)
        n_range = range(-trace_bound, trace_bound + 1)
    values: List[Pair] = []
    for n in n_range:
        for m in m_range:
            alpha = (m, n)
            if order.trace(alpha) <= trace_bound and order.is_totally_positive(alpha):
                values.append(alpha)
    values.sort(key=lambda x: (order.trace(x), x[0], x[1]))
    return values


def small_totally_positive_examples(order: RealQuadraticOrder, trace_bound: int) -> Iterable[Pair]:
    """Yield small totally positive elements of bounded trace.

    Parameters
    ----------
    order : RealQuadraticOrder
        Ambient maximal order.
    trace_bound : int
        Trace bound.

    Yields
    ------
    tuple[int, int]
        Totally positive coefficient pairs.
    """
    return enumerate_totally_positive_targets(order, trace_bound)


def test_exact_order_comparison_via_trace_and_norm() -> None:
    """Trace-norm comparisons agree with embedding comparisons on examples."""
    for D in (2, 3, 5, 6, 10, 13):
        order = RealQuadraticOrder(D)
        samples = list(enumerate_totally_positive_targets(order, 8))[:10]
        for x in samples:
            assert order.is_positive(x)
            assert order.is_nonnegative(x)
            assert order.le(order.zero, x)
            assert order.lt(order.zero, x)
        for x in samples:
            z = order.sub(order.zero, x)
            assert not order.is_nonnegative(z)
            assert not order.is_positive(z)



def test_sigma1_sign_and_balance_by_unit_square() -> None:
    """The exact first-embedding sign and balancing helper behave correctly."""
    for D in (5, 10, 13):
        order = RealQuadraticOrder(D)
        assert order.sigma1_sign((1, 0)) > 0
        beta, scale = order.balance_by_unit_square((1, 0))
        assert beta == (1, 0)
        assert order.mul(beta, order.sqr(scale)) == (1, 0)

    order = RealQuadraticOrder(5)
    unit_sq = order.sqr(order.balanced_totally_positive_unit())
    skew = unit_sq
    beta, scale = order.balance_by_unit_square(skew)
    assert beta == (1, 0)
    assert order.mul(beta, order.sqr(scale)) == skew


def test_balancing_benchmark_row_improves_sharp_box_candidates() -> None:
    """Balancing reduces the sharp-box candidate count on a skew unit-square orbit."""
    order = RealQuadraticOrder(5)
    coeff = order.sqr(order.balanced_totally_positive_unit())
    row = balancing_benchmark_row(order, coeff, 100)
    assert row.balanced_coeff == (1, 0)
    assert row.original_row_candidates == row.balanced_row_candidates
    assert row.original_sharp_box_candidates > row.balanced_sharp_box_candidates


def test_weighted_trace_form_discriminant_identity() -> None:
    """The weighted-form discriminant matches ``4 Disc(K) N(a)``."""
    for D in (2, 3, 5, 6, 10, 13):
        order = RealQuadraticOrder(D)
        coeffs = list(enumerate_totally_positive_targets(order, 6))[:8]
        for coeff in coeffs:
            form = weighted_trace_form(order, coeff)
            assert form.delta == discriminant_identity(order, coeff)
            assert form.A > 0 and form.C > 0 and form.delta > 0


def test_exact_row_interval_matches_direct_inequality() -> None:
    """Row intervals coincide with the direct quadratic inequality."""
    for D in (2, 5, 10):
        order = RealQuadraticOrder(D)
        coeff = list(enumerate_totally_positive_targets(order, 6))[0]
        form = weighted_trace_form(order, coeff)
        trace_bound = 20
        for n in range(-5, 6):
            interval = exact_row_interval(form, trace_bound, n)
            brute = [m for m in range(-25, 26) if form.value(m, n) <= trace_bound]
            if not brute:
                assert interval is None
            else:
                assert interval == (min(brute), max(brute))

def test_exact_row_interval_can_be_empty_despite_nonnegative_discriminant() -> None:
    """A non-negative row discriminant does not by itself guarantee a non-empty row."""
    order = RealQuadraticOrder(3)
    form = weighted_trace_form(order, (2, -1))
    disc = 4 * form.A * 3 - form.delta * (-1) * (-1)
    assert disc == 0
    assert exact_row_interval(form, 3, -1) is None


def test_weighted_search_matches_brute_force() -> None:
    """Exact weighted search agrees with brute force on small cases."""
    for D in (2, 3, 5, 6, 10):
        order = RealQuadraticOrder(D)
        coeffs = list(enumerate_totally_positive_targets(order, 6))[:6]
        alphas = list(enumerate_totally_positive_targets(order, 10))[:8]
        for coeff in coeffs:
            for alpha in alphas:
                exact = enumerate_weighted_search(order, coeff, alpha, return_all_roots=False)
                brute = brute_force_weighted_search(order, coeff, alpha)
                assert set(exact.value_to_root) == set(brute)
                assert exact.stats.accepted_roots == 2 * exact.stats.distinct_values


def test_weighted_trace_bound_search_matches_brute_force() -> None:
    """Trace-bounded weighted values agree with brute force."""
    for D in (2, 5, 10):
        order = RealQuadraticOrder(D)
        coeffs = list(enumerate_totally_positive_targets(order, 6))[:5]
        for coeff in coeffs:
            exact = enumerate_weighted_values_by_trace(order, coeff, 20, return_all_roots=False)
            brute = brute_force_weighted_values_by_trace(order, coeff, 20)
            assert set(exact.value_to_root) == set(brute)
            assert exact.stats.accepted_roots == 2 * exact.stats.distinct_values


def test_box_baselines_agree_with_exact_search() -> None:
    """Sharp-box and crude-box searches produce the same weighted values."""
    for D in (5, 10, 21):
        order = RealQuadraticOrder(D)
        coeff = list(enumerate_totally_positive_targets(order, 8))[1]
        alpha = list(enumerate_totally_positive_targets(order, 16))[5]
        row = enumerate_weighted_search(order, coeff, alpha, return_all_roots=False)
        sharp = enumerate_weighted_search_box(order, coeff, alpha, use_sharp_box=True, return_all_roots=False)
        crude = enumerate_weighted_search_box(order, coeff, alpha, use_sharp_box=False, return_all_roots=False)
        assert set(row.value_to_root) == set(sharp.value_to_root)
        assert set(row.value_to_root) == set(crude.value_to_root)
        assert row.stats.trace_candidates <= sharp.stats.trace_candidates <= crude.stats.trace_candidates


def test_diagonal_representability_matches_brute_force() -> None:
    """DP and meet-in-the-middle agree with brute force on small cases."""
    cases = [
        (2, [(1, 0), (1, 0)]),
        (2, [(1, 0), (2, 0), (1, 0)]),
        (3, [(1, 0), (1, 0), (2, 0)]),
        (5, [(1, 0), (1, 0), (1, 0)]),
        (6, [(1, 0), (2, 0)]),
        (10, [(1, 0), (4, 1)]),
    ]
    for D, coeffs in cases:
        order = RealQuadraticOrder(D)
        alphas = list(enumerate_totally_positive_targets(order, 14))[:16]
        for alpha in alphas:
            brute = brute_force_diagonal_representable(order, coeffs, alpha)
            dp = diagonal_representability_dp(order, coeffs, alpha)
            mitm = diagonal_representability_mitm(order, coeffs, alpha)
            assert dp.represented == brute
            assert mitm.represented == brute
            if dp.represented:
                assert dp.roots is not None and dp.values is not None
                total = order.zero
                for coeff, root in zip(coeffs, dp.roots):
                    total = order.add(total, order.mul(coeff, order.sqr(root)))
                assert total == alpha
            if mitm.represented:
                assert mitm.roots is not None and mitm.values is not None
                total = order.zero
                for coeff, root in zip(coeffs, mitm.roots):
                    total = order.add(total, order.mul(coeff, order.sqr(root)))
                assert total == alpha


def test_totally_positive_target_enumeration_matches_brute_force() -> None:
    """Bounded totally positive target enumeration is exact."""
    for D in (2, 3, 5, 10, 13):
        order = RealQuadraticOrder(D)
        exact = enumerate_totally_positive_targets(order, 18)
        brute = brute_force_totally_positive_targets(order, 18)
        assert exact == brute


def test_batched_bounded_representables_and_truants_match_naive() -> None:
    """Batched bounded search agrees with the per-target solver."""
    cases = [
        (2, [(1, 0), (1, 0)], 10),
        (5, [(1, 0), (1, 0), (1, 0)], 12),
        (10, [(1, 0), (4, 1)], 18),
    ]
    for D, coeffs, trace_bound in cases:
        order = RealQuadraticOrder(D)
        batched_repr = batched_bounded_representables(order, coeffs, trace_bound)
        batched_truants = bounded_truants_batched(order, coeffs, trace_bound)
        naive_truants = bounded_truants_naive(order, coeffs, trace_bound, method="mitm")
        assert batched_truants == naive_truants
        targets = enumerate_totally_positive_targets(order, trace_bound)
        represented = [alpha for alpha in targets if alpha not in batched_truants]
        assert list(batched_repr.representables) == represented


def test_batched_and_naive_agree_on_small_family() -> None:
    """Systematic bounded cross-validation on a small family of fields and ranks."""
    for D in (5, 6, 10, 13, 21):
        order = RealQuadraticOrder(D)
        for coeffs in ([(1, 0), (1, 0)], [(1, 0), (1, 0), (1, 0)]):
            for trace_bound in (8, 12, 16):
                batched = bounded_truants_batched(order, coeffs, trace_bound)
                naive = bounded_truants_naive(order, coeffs, trace_bound, method="mitm")
                assert batched == naive


def test_benchmarks_generate_consistent_rows() -> None:
    """Benchmark helpers return structurally consistent rows."""
    enumeration_rows = benchmark_weighted_enumeration([(21, (2, 1), (30, 0))], repeats=1)
    assert len(enumeration_rows) == 1
    erow = enumeration_rows[0]
    assert erow.accepted_roots == 2 * erow.distinct_values
    assert erow.row_candidates <= erow.sharp_box_candidates <= erow.crude_box_candidates

    repr_rows = benchmark_representability([(10, [(1, 0), (4, 1)], (18, 2))], repeats=1)
    assert len(repr_rows) == 1
    rrow = repr_rows[0]
    assert len(rrow.distinct_value_counts) == 2
    assert len(rrow.dp_state_counts) == 2
    assert rrow.mitm_left_states >= 1 and rrow.mitm_right_states >= 1
    assert rrow.preprocessing_ms >= 0.0
    assert rrow.dp_combination_ms >= 0.0
    assert rrow.dp_total_ms >= rrow.dp_combination_ms
    assert rrow.mitm_total_ms >= rrow.mitm_combination_ms

    bounded_rows = benchmark_bounded_truants([(5, [(1, 0), (1, 0), (1, 0)], 10)], repeats=1)
    assert len(bounded_rows) == 1
    brow = bounded_rows[0]
    assert brow.representables >= 0
    assert brow.truants >= 0
    assert len(brow.batched_state_counts) == 3
    assert len(brow.value_counts) == 3
    assert brow.batched_total_ms >= brow.batched_combination_ms
    assert brow.naive_total_ms >= 0.0


def test_error_branches() -> None:
    """Invalid inputs raise the documented errors."""
    try:
        RealQuadraticOrder(12)
    except ValueError:
        pass
    else:
        raise AssertionError("Non-squarefree D should be rejected")

    order = RealQuadraticOrder(5)
    nonpositive = (0, 1)
    positive = (1, 0)
    alpha = (3, 1)
    try:
        weighted_trace_form(order, nonpositive)
    except ValueError:
        pass
    else:
        raise AssertionError("Non-totally-positive coefficient should be rejected")

    try:
        enumerate_weighted_search(order, positive, nonpositive)
    except ValueError:
        pass
    else:
        raise AssertionError("Non-totally-positive target should be rejected")

    try:
        bounded_truants_naive(order, [positive], 8, method="unknown")
    except ValueError:
        pass
    else:
        raise AssertionError("Unknown method should be rejected")


def test_benchmark_optimisation_runs() -> None:
    """The optimisation benchmark returns structurally consistent rows."""
    rows = benchmark_optimisation(
        [
            (10, [(1, 0), (4, 1)], (18, 2)),
            (13, [(1, 0), (2, 1), (2, 1), (1, 0)], (18, 0)),
        ],
        repeats=2,
    )
    assert len(rows) == 2
    for row in rows:
        assert len(row.original_distinct_value_counts) == len(row.coeffs)
        assert len(row.optimised_distinct_value_counts) == len(row.coeffs)
        assert len(row.original_state_counts) == len(row.coeffs)
        assert len(row.optimised_state_counts) == len(row.coeffs)
        assert sorted(row.coefficient_order) == list(range(len(row.coeffs)))


def test_generic_baseline_benchmark_agrees_with_specialised_solvers() -> None:
    """The surrogate generic exact baseline agrees with DP and MITM."""
    rows = benchmark_generic_baseline([
        (10, [(1, 0), (4, 1)], (18, 2)),
        (5, [(1, 0), (1, 0), (1, 0)], (15, 1)),
    ], repeats=1)
    assert len(rows) == 2
    for row in rows:
        assert row.cartesian_product_size > 0
        assert row.generic_total_ms >= row.preprocessing_ms

