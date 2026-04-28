#!/usr/bin/env python3
"""Regenerate table data and performance plots for the quadratic_sos paper."""

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


def _format_counts(D: int, trace_bound: int) -> str:
    order = RealQuadraticOrder(D)
    counts, first = length_distribution(order, trace_bound)
    keys = [1, 2, 3, 4, 5, None]
    lines = [f"D={D}", f"trace_bound={trace_bound}"]
    for key in keys:
        lines.append(f"length_{key if key is not None else 'infty'}={counts.get(key, 0)}")
        if key in first:
            lines.append(f"first_{key if key is not None else 'infty'}={first[key]}")
    return "\n".join(lines) + "\n"


def _cross_field(fields: list[int], trace_bound: int) -> str:
    lines = ["D\ttotal\tlen1\tlen2\tlen3\tlen4\tlen5\tinfty"]
    for D in fields:
        order = RealQuadraticOrder(D)
        counts, _ = length_distribution(order, trace_bound)
        total = sum(counts.values())
        row = [D, total, counts.get(1, 0), counts.get(2, 0), counts.get(3, 0), counts.get(4, 0), counts.get(5, 0), counts.get(None, 0)]
        lines.append("\t".join(map(str, row)))
    return "\n".join(lines) + "\n"


def _balancing(D: int, beta: tuple[int, int], max_k: int) -> tuple[str, list[tuple[int, int, int]]]:
    order = RealQuadraticOrder(D)
    lines = ["k\traw_trace\traw_trace_v_bound\traw_rows\traw_box\tbalanced_trace\tbalanced_trace_v_bound\tbalanced_rows\tbalanced_box"]
    plot_rows: list[tuple[int, int, int]] = []
    for k, alpha in enumerate(pell_orbit(order, beta, max_k)):
        raw = search_profile(order, alpha, balance=False)
        bal = search_profile(order, alpha, balance=True)
        lines.append(
            f"{k}\t{raw.trace}\t{raw.trace_v_bound}\t{raw.scanned_rows}\t{raw.box_size}\t"
            f"{bal.trace}\t{bal.trace_v_bound}\t{bal.scanned_rows}\t{bal.box_size}"
        )
        plot_rows.append((k, raw.trace_v_bound, bal.trace_v_bound))
    return "\n".join(lines) + "\n", plot_rows


def _plot_balancing(plot_rows: list[tuple[int, int, int]]) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - optional dependency
        print(f"matplotlib unavailable, skipping plot: {exc}")
        return
    fig_dir = ROOT / "paper" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    xs = [row[0] for row in plot_rows]
    raw = [row[1] for row in plot_rows]
    bal = [row[2] for row in plot_rows]
    plt.figure(figsize=(6.0, 3.4))
    plt.plot(xs, raw, marker="o", linestyle="-", color="black", label="raw")
    plt.plot(xs, bal, marker="s", linestyle="--", color="black", label="balanced")
    plt.xlabel("Pell-orbit index k")
    plt.ylabel("trace-row bound")
    plt.yscale("log")
    plt.legend()
    plt.tight_layout()
    plt.savefig(fig_dir / "sos_balancing_rows.pdf")
    plt.savefig(fig_dir / "sos_balancing_rows.png", dpi=180)
    plt.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true", help="use smaller bounds for fast smoke testing")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data = ROOT / "data"
    data.mkdir(exist_ok=True)
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
    (data / "sos_D10_distribution.txt").write_text(_format_counts(10, trace_bound))
    (data / "sos_cross_field.tsv").write_text(_cross_field(fields, cross_bound))
    balancing_text, plot_rows = _balancing(10, (18, 2), max_k)
    (data / "sos_balancing_ablation.tsv").write_text(balancing_text)
    _plot_balancing(plot_rows)
    summary = validation_sweep(val_max_D, val_trace)
    (data / "sos_validation_summary.tsv").write_text(
        "fields\ttrace_bound\telements\tmismatches\n"
        f"{summary.field_count}\t{summary.trace_bound}\t{summary.element_count}\t{summary.mismatches}\n"
    )
    print(f"Wrote data to {data}")
    print(f"Validation mismatches: {summary.mismatches}")


if __name__ == "__main__":
    main()
