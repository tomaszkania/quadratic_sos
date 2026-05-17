# TOMS Algorithm-paper software component

Version: **0.2.0**.

This repository is intended to be submitted as the software component of a TOMS Algorithm paper.
The paper PDF is submitted separately; the `paper/` directory also contains the ACM `acmart`
LaTeX source, BibTeX file, generated figures, and a locally compiled PDF for convenience.
The directory `submission/` contains a generated convenience archive with the same software
component.

## Required commands

The following commands require no interactive input.

```bash
python -m pip install -e .[repro]
python scripts/run_all_checks.py
```

The smoke-test driver runs the regression suite, selected reproduction scripts, and validation
checks. Timing columns are treated as machine-dependent; deterministic structural counts are
regenerated from the package code.

## Contents

- `src/`: installable Python package.
- `tests/`: regression and brute-force cross-checks.
- `scripts/run_all_checks.py`: non-interactive smoke-test driver.
- `scripts/reproduce_tables.py`: regenerates table data and figures.
- `scripts/validation_sweep.py`: bounded validation driver.
- `scripts/make_submission_archive.py`: builds the TOMS software-component archive.
- `data/`: generated deterministic data and timing summaries.
- `notebooks/`: executed paper-illustration notebook.
- `paper/`: ACM manuscript source, BibTeX references, generated figures, and PDF.

## Portability

The core package uses only the Python standard library and arbitrary-precision integers. The
optional `repro` extra installs `pytest`, `matplotlib`, `jupyter`, and `nbconvert` for tests,
figures, and notebooks. The package metadata declares Python 3.10 or later; this bundle was
smoke-tested with Python 3.13.

## Expected outputs

A successful smoke test ends by writing `data/run_all_checks_summary.txt` and printing a final
`OK` line. The number of regression tests is recorded in that summary file.
