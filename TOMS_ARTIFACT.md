# TOMS software component

This repository is the software component accompanying the TOMS submission
*Computing Integral Length of Sums of Squares and Pythagoras Elements in Real Quadratic Rings*.
It is intended for evaluation as a documented, complete, and exercisable mathematical-software
artefact.

## Required commands

The following commands require no interactive input:

```bash
python -m pip install -e .[repro]
python scripts/run_all_checks.py
```

A successful smoke test writes `data/run_all_checks_summary.txt` and prints a final `OK` line.
Quick smoke-test data are written under `data/smoke/` so that the deterministic paper data in
`data/` are not overwritten.

## Full reproduction commands

```bash
python scripts/reproduce_tables.py
python scripts/validation_sweep.py --max-D 41 --trace-bound 24
```

Timing values in the balancing-ablation script are machine-dependent. Structural counts, exact
lengths, validation counts, and witness decompositions are deterministic.

## Contents

- `src/`: installable Python package using only standard-library runtime dependencies.
- `tests/`: regression tests for arithmetic, representability, exact lengths, decompositions,
  Pythagoras witnesses, balancing, validation, and table counts.
- `scripts/run_all_checks.py`: non-interactive smoke-test driver.
- `scripts/reproduce_tables.py`: deterministic table-data and figure-generation driver.
- `scripts/validation_sweep.py`: bounded brute-force validation driver.
- `data/`: deterministic data corresponding to the computational section.
- `notebooks/`: executed illustration notebook.
- `paper/`: ACM/TOMS LaTeX source and compiled PDF.

## Portability

The core package uses Python arbitrary-precision integers and no external runtime library. The
optional `repro` extra installs `pytest`, `matplotlib`, `jupyter`, and `nbconvert` for tests,
figures, and notebook execution. The smoke-test path has been checked in this bundle with
Python 3.13.

## Validation independence

The validation sweep compares the main exact-length routine against a deliberately naive
coefficient-box search. That validator does not call Peters' representability test, does not use the
row-wise search engine, and does not use the closed-form Pythagoras-number classification; it tests
all sums of at most five relevant squares inside the coefficient box.
