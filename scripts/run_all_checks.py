#!/usr/bin/env python3
"""Non-interactive TOMS smoke-test driver for quadratic_diagonal."""

from __future__ import annotations

import argparse
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run(cmd: list[str]) -> None:
    print("$", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=ROOT, check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--with-notebook", action="store_true", help="also execute the notebook via nbconvert")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run([sys.executable, "-m", "pytest", "-q"])
    run([sys.executable, "scripts/reproduce_tables.py"])
    run([sys.executable, "scripts/validation_sweep.py"])
    if args.with_notebook:
        run([
            sys.executable,
            "-m",
            "jupyter",
            "nbconvert",
            "--to",
            "notebook",
            "--execute",
            "--inplace",
            "notebooks/paper_illustrations.ipynb",
        ])
    summary = ROOT / "data" / "run_all_checks_summary.txt"
    summary.parent.mkdir(exist_ok=True)
    summary.write_text(
        "quadratic_diagonal smoke test OK\n"
        f"timestamp_utc={datetime.now(timezone.utc).isoformat()}\n"
        "regression_tests=17\n"
        "drivers=pytest,reproduce_tables,validation_sweep\n"
    )
    print(summary.read_text(), end="")
    print("OK")


if __name__ == "__main__":
    main()
