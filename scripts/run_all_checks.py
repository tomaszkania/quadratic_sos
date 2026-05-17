#!/usr/bin/env python3
"""Non-interactive TOMS smoke-test driver for quadratic_sos."""

from __future__ import annotations

import argparse
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from quadratic_sos import __version__


def run(cmd: list[str]) -> None:
    """Run a repository command and fail immediately on nonzero exit."""
    print("$", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=ROOT, check=True)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--with-notebook", action="store_true", help="also execute the notebook via nbconvert")
    return parser.parse_args()


def main() -> None:
    """Run the smoke-test workflow."""
    args = parse_args()
    run([sys.executable, "-m", "pytest", "-q"])
    run([sys.executable, "scripts/reproduce_tables.py", "--quick", "--output-dir", "data/smoke"])
    run([sys.executable, "scripts/validation_sweep.py", "--max-D", "17", "--trace-bound", "16"])
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
        "quadratic_sos smoke test OK\n"
        f"version={__version__}\n"
        f"timestamp_utc={datetime.now(timezone.utc).isoformat()}\n"
        "regression_tests=26\n"
        "drivers=pytest,reproduce_tables_quick,validation_sweep_small\n",
        encoding="utf-8",
    )
    print(summary.read_text(encoding="utf-8"), end="")
    print("OK")


if __name__ == "__main__":
    main()
