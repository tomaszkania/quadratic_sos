"""Generate the paper tables and performance plots.

This script regenerates the TSV tables in ``data/`` and the performance plots in
``paper/figures/`` from the public package API.
"""

from __future__ import annotations

from pathlib import Path
import sys
from typing import Any, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from quadratic_diagonal import (
    RealQuadraticOrder,
    balancing_benchmark_row,
    benchmark_bounded_truants,
    benchmark_generic_baseline,
    benchmark_optimisation,
    benchmark_representability,
    benchmark_weighted_enumeration,
    bounded_truants_batched,
    format_pair,
)


def _write_enumeration_table(path: Path) -> list[Any]:
    """Write the weighted-enumeration benchmark table.

    Parameters
    ----------
    path : pathlib.Path
        Output TSV path.

    Returns
    -------
    list
        Benchmark rows as returned by :func:`benchmark_weighted_enumeration`.
    """
    rows = benchmark_weighted_enumeration(
        [
            (21, (2, 1), (30, 0)),
            (21, (2, 1), (40, 0)),
            (6, (3, -1), (30, 0)),
            (5, (1, 1), (30, 0)),
            (10, (1, 0), (40, 0)),
        ],
        repeats=5,
    )
    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            "D\tcoeff\talpha\tdistinct_values\ttotal_roots\tsharp_box_candidates\t"
            "crude_box_candidates\trow_candidates\tsharp_box_ms\tcrude_box_ms\trow_ms\n"
        )
        for row in rows:
            handle.write(
                f"{row.D}\t{format_pair(row.coeff)}\t{format_pair(row.alpha)}\t{row.distinct_values}\t"
                f"{row.accepted_roots}\t{row.sharp_box_candidates}\t{row.crude_box_candidates}\t"
                f"{row.row_candidates}\t{row.sharp_box_ms:.3f}\t{row.crude_box_ms:.3f}\t{row.row_ms:.3f}\n"
            )
    return rows


def _write_representability_table(path: Path) -> list[Any]:
    """Write the representability benchmark table.

    Parameters
    ----------
    path : pathlib.Path
        Output TSV path.

    Returns
    -------
    list
        Benchmark rows as returned by :func:`benchmark_representability`.
    """
    rows = benchmark_representability(
        [
            (10, [(1, 0), (4, 1)], (18, 2)),
            (5, [(1, 0), (1, 0), (1, 0)], (15, 1)),
            (21, [(2, 1), (1, 0), (1, 0), (1, 0)], (30, 0)),
            (13, [(1, 0), (2, 1), (2, 1), (1, 0)], (18, 0)),
        ],
        repeats=5,
    )
    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            "D\tcoeffs\talpha\tinternal_order\tdistinct_value_counts\tdp_state_counts\t"
            "mitm_left_states\tmitm_right_states\tpreprocessing_ms\tdp_combination_ms\t"
            "dp_total_ms\tmitm_combination_ms\tmitm_total_ms\n"
        )
        for row in rows:
            coeffs = "[" + ",".join(format_pair(x) for x in row.coeffs) + "]"
            values = "[" + ",".join(str(x) for x in row.distinct_value_counts) + "]"
            dp_states = "[" + ",".join(str(x) for x in row.dp_state_counts) + "]"
            internal_order = "[" + ",".join(str(x) for x in getattr(row, "coefficient_order", ())) + "]"
            handle.write(
                f"{row.D}\t{coeffs}\t{format_pair(row.alpha)}\t{internal_order}\t{values}\t{dp_states}\t"
                f"{row.mitm_left_states}\t{row.mitm_right_states}\t{row.preprocessing_ms:.3f}\t"
                f"{row.dp_combination_ms:.3f}\t{row.dp_total_ms:.3f}\t"
                f"{row.mitm_combination_ms:.3f}\t{row.mitm_total_ms:.3f}\n"
            )
    return rows


def _write_generic_baseline_table(path: Path) -> list[Any]:
    """Write the surrogate general-purpose baseline benchmark table.

    Parameters
    ----------
    path : pathlib.Path
        Output TSV path.

    Returns
    -------
    list
        Benchmark rows as returned by :func:`benchmark_generic_baseline`.
    """
    rows = benchmark_generic_baseline(
        [
            (10, [(1, 0), (4, 1)], (18, 2)),
            (5, [(1, 0), (1, 0), (1, 0)], (15, 1)),
            (21, [(2, 1), (1, 0), (1, 0), (1, 0)], (30, 0)),
            (13, [(1, 0), (2, 1), (2, 1), (1, 0)], (18, 0)),
        ],
        repeats=3,
    )
    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            "D\tcoeffs\talpha\tinternal_order\tdistinct_value_counts\tcartesian_product_size\t"
            "preprocessing_ms\tdp_total_ms\tmitm_total_ms\tgeneric_total_ms\n"
        )
        for row in rows:
            coeffs = "[" + ",".join(format_pair(x) for x in row.coeffs) + "]"
            internal_order = "[" + ",".join(str(x) for x in row.coefficient_order) + "]"
            value_counts = "[" + ",".join(str(x) for x in row.distinct_value_counts) + "]"
            handle.write(
                f"{row.D}\t{coeffs}\t{format_pair(row.alpha)}\t{internal_order}\t{value_counts}\t"
                f"{row.cartesian_product_size}\t{row.preprocessing_ms:.3f}\t{row.dp_total_ms:.3f}\t"
                f"{row.mitm_total_ms:.3f}\t{row.generic_total_ms:.3f}\n"
            )
    return rows


def _write_bounded_table(path: Path) -> list[Any]:
    """Write the bounded-search benchmark table.

    Parameters
    ----------
    path : pathlib.Path
        Output TSV path.

    Returns
    -------
    list
        Benchmark rows as returned by :func:`benchmark_bounded_truants`.
    """
    rows = benchmark_bounded_truants(
        [
            (5, [(1, 0), (1, 0), (1, 0)], 10),
            (10, [(1, 0), (4, 1)], 18),
            (21, [(2, 1), (1, 0), (1, 0)], 16),
        ],
        repeats=3,
    )
    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            "D\tcoeffs\ttrace_bound\trepresentables\ttruants\tbatched_state_counts\tvalue_counts\t"
            "batched_preprocessing_ms\tbatched_combination_ms\tbatched_total_ms\tnaive_total_ms\n"
        )
        for row in rows:
            coeffs = "[" + ",".join(format_pair(x) for x in row.coeffs) + "]"
            states = "[" + ",".join(str(x) for x in row.batched_state_counts) + "]"
            value_counts = "[" + ",".join(str(x) for x in row.value_counts) + "]"
            handle.write(
                f"{row.D}\t{coeffs}\t{row.trace_bound}\t{row.representables}\t{row.truants}\t"
                f"{states}\t{value_counts}\t{row.batched_preprocessing_ms:.3f}\t"
                f"{row.batched_combination_ms:.3f}\t{row.batched_total_ms:.3f}\t{row.naive_total_ms:.3f}\n"
            )
    return rows


def _write_optimisation_table(path: Path) -> list[Any]:
    """Write the optimisation benchmark table.

    Parameters
    ----------
    path : pathlib.Path
        Output TSV path.

    Returns
    -------
    list
        Benchmark rows as returned by :func:`benchmark_optimisation`.
    """
    rows = benchmark_optimisation(
        [
            (10, [(1, 0), (4, 1)], (18, 2)),
            (13, [(1, 0), (2, 1), (2, 1), (1, 0)], (18, 0)),
            (21, [(2, 1), (1, 0), (1, 0), (1, 0)], (30, 0)),
        ],
        repeats=5,
    )
    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            "D\tcoeffs\talpha\tinternal_order\toriginal_value_counts\toptimised_value_counts\t"
            "original_state_counts\toptimised_state_counts\tuncached_preprocessing_ms\t"
            "cached_preprocessing_ms\toriginal_total_ms\toptimised_total_ms\n"
        )
        for row in rows:
            coeffs = "[" + ",".join(format_pair(x) for x in row.coeffs) + "]"
            internal_order = "[" + ",".join(str(x) for x in row.coefficient_order) + "]"
            original_values = "[" + ",".join(str(x) for x in row.original_distinct_value_counts) + "]"
            optimised_values = "[" + ",".join(str(x) for x in row.optimised_distinct_value_counts) + "]"
            original_states = "[" + ",".join(str(x) for x in row.original_state_counts) + "]"
            optimised_states = "[" + ",".join(str(x) for x in row.optimised_state_counts) + "]"
            handle.write(
                f"{row.D}\t{coeffs}\t{format_pair(row.alpha)}\t{internal_order}\t{original_values}\t"
                f"{optimised_values}\t{original_states}\t{optimised_states}\t"
                f"{row.uncached_preprocessing_ms:.3f}\t{row.cached_preprocessing_ms:.3f}\t"
                f"{row.original_total_ms:.3f}\t{row.optimised_total_ms:.3f}\n"
            )
    return rows


def _write_balancing_table(path: Path) -> list[Any]:
    """Write the balancing benchmark table.

    Parameters
    ----------
    path : pathlib.Path
        Output TSV path.

    Returns
    -------
    list
        Singleton list containing the balancing benchmark row.
    """
    order = RealQuadraticOrder(5)
    coeff = order.sqr(order.balanced_totally_positive_unit())
    row = balancing_benchmark_row(order, coeff, 100)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            "D\tcoeff\tbalanced_coeff\ttrace_bound\toriginal_row_candidates\tbalanced_row_candidates\t"
            "original_sharp_box_candidates\tbalanced_sharp_box_candidates\n"
        )
        handle.write(
            f"{row.D}\t{format_pair(row.coeff)}\t{format_pair(row.balanced_coeff)}\t{row.trace_bound}\t"
            f"{row.original_row_candidates}\t{row.balanced_row_candidates}\t"
            f"{row.original_sharp_box_candidates}\t{row.balanced_sharp_box_candidates}\n"
        )
    return [row]


def _write_arithmetic_examples(path: Path) -> None:
    """Write the arithmetic examples table.

    Parameters
    ----------
    path : pathlib.Path
        Output TSV path.
    """
    cases = [
        (5, [(1, 0), (1, 0), (1, 0)], 12),
        (5, [(1, 0), (1, 0), (2, 1)], 40),
        (13, [(1, 0), (2, 1), (2, 1), (1, 0)], 18),
        (21, [(2, 1), (1, 0), (1, 0)], 16),
    ]
    with path.open("w", encoding="utf-8") as handle:
        handle.write("D\tcoeffs\ttrace_bound\ttruants\n")
        for D, coeffs, B in cases:
            order = RealQuadraticOrder(D)
            truants = bounded_truants_batched(order, coeffs, B)
            coeffs_text = "[" + ",".join(format_pair(x) for x in coeffs) + "]"
            truants_text = "[" + ",".join(format_pair(x) for x in truants) + "]"
            handle.write(f"{D}\t{coeffs_text}\t{B}\t{truants_text}\n")


def _hatched_bar(
    ax: Any,
    positions: Sequence[float],
    values: Sequence[float],
    *,
    width: float,
    label: str,
    hatch: str,
) -> Any:
    """Draw a monochrome-friendly hatched bar group."""
    bars = ax.bar(positions, values, width=width, label=label, facecolor="white", edgecolor="black")
    for bar in bars:
        bar.set_hatch(hatch)
    return bars


def _save_plot(fig: plt.Figure, pdf_path: Path, png_path: Path) -> None:
    """Save a Matplotlib figure in both PDF and PNG formats.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to save.
    pdf_path : pathlib.Path
        PDF output path.
    png_path : pathlib.Path
        PNG output path.
    """
    fig.tight_layout()
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_enumeration_candidates(rows: Sequence[Any], out_dir: Path) -> None:
    """Plot canonical candidate counts for enumeration benchmarks.

    Parameters
    ----------
    rows : sequence
        Enumeration benchmark rows.
    out_dir : pathlib.Path
        Output figure directory.
    """
    labels = [f"D={row.D}\n{format_pair(row.coeff)}\n{format_pair(row.alpha)}" for row in rows]
    x = list(range(len(rows)))
    width = 0.25
    fig = plt.figure(figsize=(8.0, 3.8))
    ax = fig.add_subplot(1, 1, 1)
    _hatched_bar(ax, [i - width for i in x], [row.row_candidates for row in rows], width=width, label="row", hatch="")
    _hatched_bar(ax, x, [row.sharp_box_candidates for row in rows], width=width, label="sharp", hatch="//")
    _hatched_bar(ax, [i + width for i in x], [row.crude_box_candidates for row in rows], width=width, label="crude", hatch="xx")
    ax.set_yscale("log")
    ax.set_ylabel("canonical candidates")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend(frameon=False, ncols=3)
    ax.set_title("Enumeration candidates by method")
    _save_plot(fig, out_dir / "enumeration_candidates.pdf", out_dir / "enumeration_candidates.png")


def _plot_generic_baseline_runtime(rows: Sequence[Any], out_dir: Path) -> None:
    """Plot end-to-end runtimes against the generic exact baseline.

    Parameters
    ----------
    rows : sequence
        Generic-baseline benchmark rows.
    out_dir : pathlib.Path
        Output figure directory.
    """
    labels = [f"D={row.D}\n{format_pair(row.alpha)}" for row in rows]
    x = list(range(len(rows)))
    width = 0.25
    fig = plt.figure(figsize=(7.0, 3.9))
    ax = fig.add_subplot(1, 1, 1)
    _hatched_bar(ax, [i - width for i in x], [row.dp_total_ms for row in rows], width=width, label="DP", hatch="")
    _hatched_bar(ax, x, [row.mitm_total_ms for row in rows], width=width, label="MITM", hatch="//")
    _hatched_bar(ax, [i + width for i in x], [row.generic_total_ms for row in rows], width=width, label="generic", hatch="xx")
    ax.set_ylabel("end-to-end runtime (ms)")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend(frameon=False, ncols=3)
    ax.set_title("Specialised solvers versus a generic exact baseline")
    _save_plot(fig, out_dir / "generic_baseline_runtime.pdf", out_dir / "generic_baseline_runtime.png")


def _plot_bounded_runtime(rows: Sequence[Any], out_dir: Path) -> None:
    """Plot end-to-end runtimes for batched versus naive bounded search.

    Parameters
    ----------
    rows : sequence
        Bounded-search benchmark rows.
    out_dir : pathlib.Path
        Output figure directory.
    """
    labels = [f"D={row.D}\nB={row.trace_bound}" for row in rows]
    x = list(range(len(rows)))
    width = 0.35
    fig = plt.figure(figsize=(6.0, 3.8))
    ax = fig.add_subplot(1, 1, 1)
    _hatched_bar(ax, [i - width / 2 for i in x], [row.batched_total_ms for row in rows], width=width, label="batched", hatch="")
    _hatched_bar(ax, [i + width / 2 for i in x], [row.naive_total_ms for row in rows], width=width, label="naive", hatch="//")
    ax.set_ylabel("end-to-end runtime (ms)")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend(frameon=False)
    ax.set_title("Batched bounded search versus naive recomputation")
    _save_plot(fig, out_dir / "bounded_runtime.pdf", out_dir / "bounded_runtime.png")


def _plot_optimisation_runtime(rows: Sequence[Any], out_dir: Path) -> None:
    """Plot end-to-end runtime before and after ordering/caching optimisation.

    Parameters
    ----------
    rows : sequence
        Optimisation benchmark rows.
    out_dir : pathlib.Path
        Output figure directory.
    """
    labels = [f"D={row.D}\n{format_pair(row.alpha)}" for row in rows]
    x = list(range(len(rows)))
    width = 0.35
    fig = plt.figure(figsize=(6.5, 3.8))
    ax = fig.add_subplot(1, 1, 1)
    _hatched_bar(ax, [i - width / 2 for i in x], [row.original_total_ms for row in rows], width=width, label="original", hatch="")
    _hatched_bar(ax, [i + width / 2 for i in x], [row.optimised_total_ms for row in rows], width=width, label="optimised", hatch="//")
    ax.set_ylabel("DP end-to-end runtime (ms)")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend(frameon=False)
    ax.set_title("Optimisation effect on constructive DP")
    _save_plot(fig, out_dir / "optimisation_runtime.pdf", out_dir / "optimisation_runtime.png")


def main() -> None:
    """Generate all TSV tables and plot figures."""
    data_dir = ROOT / "data"
    fig_dir = ROOT / "paper" / "figures"
    data_dir.mkdir(exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    enum_rows = _write_enumeration_table(data_dir / "benchmark_enumeration.tsv")
    repr_rows = _write_representability_table(data_dir / "benchmark_representability.tsv")
    generic_rows = _write_generic_baseline_table(data_dir / "benchmark_generic_baseline.tsv")
    bounded_rows = _write_bounded_table(data_dir / "benchmark_bounded.tsv")
    opt_rows = _write_optimisation_table(data_dir / "benchmark_optimisation.tsv")
    _ = _write_balancing_table(data_dir / "benchmark_balancing.tsv")
    _write_arithmetic_examples(data_dir / "arithmetic_examples.tsv")

    _plot_enumeration_candidates(enum_rows, fig_dir)
    _plot_generic_baseline_runtime(generic_rows, fig_dir)
    _plot_bounded_runtime(bounded_rows, fig_dir)
    _plot_optimisation_runtime(opt_rows, fig_dir)


if __name__ == "__main__":
    main()
