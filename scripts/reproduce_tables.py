#!/usr/bin/env python3
"""Regenerate deterministic table data and the balancing plot."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from quadratic_sos import RealQuadraticOrder
from quadratic_sos.experiments import pell_orbit, search_profile, validation_sweep
from quadratic_sos.tables import length_distribution

Pair = tuple[int, int]


def _length_label(length: int | None) -> str:
    """Return a stable text label for a finite length or infinity."""
    return str(length) if length is not None else "infty"


def _format_counts(radicand: int, trace_bound: int) -> str:
    """Return the local exact-length distribution as key-value text."""
    order = RealQuadraticOrder(radicand)
    counts, first = length_distribution(order, trace_bound)
    lines = [f"D={radicand}", f"trace_bound={trace_bound}"]
    for length in [1, 2, 3, 4, 5, None]:
        label = _length_label(length)
        lines.append(f"length_{label}={counts.get(length, 0)}")
        if length in first:
            lines.append(f"first_{label}={first[length]}")
    return "\n".join(lines) + "\n"


def _cross_field(fields: list[int], trace_bound: int) -> str:
    """Return cross-field exact-length distributions as TSV text."""
    header = "D\ttotal\tlen1\tlen2\tlen3\tlen4\tlen5\tinfty\ttau_max"
    lines = [header]
    for radicand in fields:
        order = RealQuadraticOrder(radicand)
        counts, first = length_distribution(order, trace_bound)
        pythagoras = order.pythagoras_number()
        tau_max = order.trace(first[pythagoras]) if pythagoras in first else "NA"
        row = [
            radicand,
            sum(counts.values()),
            counts.get(1, 0),
            counts.get(2, 0),
            counts.get(3, 0),
            counts.get(4, 0),
            counts.get(5, 0),
            counts.get(None, 0),
            tau_max,
        ]
        lines.append("\t".join(map(str, row)))
    return "\n".join(lines) + "\n"


def _balancing(radicand: int, beta: Pair, max_k: int) -> tuple[str, list[tuple[int, int, int]]]:
    """Return the deterministic balancing-ablation data as TSV text."""
    order = RealQuadraticOrder(radicand)
    lines = [
        "k\traw_trace\traw_trace_v_bound\traw_rows\traw_box\traw_squares\t"
        "balanced_trace\tbalanced_trace_v_bound\tbalanced_rows\tbalanced_box\tbalanced_squares"
    ]
    plot_rows: list[tuple[int, int, int]] = []
    for k, alpha in enumerate(pell_orbit(order, beta, max_k)):
        raw = search_profile(order, alpha, balance=False)
        balanced = search_profile(order, alpha, balance=True)
        lines.append(
            f"{k}\t{raw.trace}\t{raw.trace_v_bound}\t{raw.scanned_rows}\t{raw.box_size}\t{raw.square_count}\t"
            f"{balanced.trace}\t{balanced.trace_v_bound}\t{balanced.scanned_rows}\t{balanced.box_size}\t{balanced.square_count}"
        )
        plot_rows.append((k, raw.trace_v_bound, balanced.trace_v_bound))
    return "\n".join(lines) + "\n", plot_rows


def _plot_balancing(plot_rows: list[tuple[int, int, int]]) -> None:
    """Regenerate the balancing-ablation plot when matplotlib is available."""
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - optional dependency
        print(f"matplotlib unavailable, skipping plot: {exc}")
        return

    fig_dir = ROOT / "paper" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    xs = [row[0] for row in plot_rows]
    raw = [row[1] for row in plot_rows]
    balanced = [row[2] for row in plot_rows]

    plt.figure(figsize=(6.0, 3.4))
    plt.plot(xs, raw, marker="o", linestyle="-", color="black", label="raw")
    plt.plot(xs, balanced, marker="s", linestyle="--", color="black", label="balanced")
    plt.xlabel("Pell-orbit index k")
    plt.ylabel("trace-row bound")
    plt.yscale("log")
    plt.legend()
    plt.tight_layout()
    plt.savefig(fig_dir / "sos_balancing_rows.pdf")
    plt.savefig(fig_dir / "sos_balancing_rows.png", dpi=180)
    plt.close()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true", help="use smaller bounds for fast smoke testing")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="directory for generated text data; defaults to data/ or data/smoke/ for --quick",
    )
    return parser.parse_args()


def main() -> None:
    """Regenerate table data."""
    args = parse_args()
    if args.quick:
        trace_bound = 24
        cross_bound = 30
        fields = [5, 10, 13]
        max_k = 2
        val_max_D = 17
        val_trace = 16
    else:
        trace_bound = 40
        cross_bound = 80
        fields = [5, 10, 13, 29, 53, 101]
        max_k = 3
        val_max_D = 41
        val_trace = 24

    if args.output_dir is not None:
        data = args.output_dir
    elif args.quick:
        data = ROOT / "data" / "smoke"
    else:
        data = ROOT / "data"
    data.mkdir(parents=True, exist_ok=True)

    (data / "sos_D10_distribution.txt").write_text(_format_counts(10, trace_bound), encoding="utf-8")
    (data / "sos_cross_field.tsv").write_text(_cross_field(fields, cross_bound), encoding="utf-8")
    balancing_text, plot_rows = _balancing(10, (18, 2), max_k)
    (data / "sos_balancing_ablation.tsv").write_text(balancing_text, encoding="utf-8")
    if not args.quick:
        _plot_balancing(plot_rows)
    summary = validation_sweep(val_max_D, val_trace)
    (data / "sos_validation_summary.tsv").write_text(
        "fields\ttrace_bound\telements\tmismatches\n"
        f"{summary.field_count}\t{summary.trace_bound}\t{summary.element_count}\t{summary.mismatches}\n",
        encoding="utf-8",
    )
    print(f"Wrote data to {data}")
    print(f"Validation mismatches: {summary.mismatches}")


if __name__ == "__main__":
    main()
