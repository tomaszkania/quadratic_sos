# TOMS submission checklist

## Article

- [x] Manuscript converted to the ACM primary article template in `manuscript` mode.
- [x] `\acmJournal{TOMS}` and ACM CCS descriptors supplied.
- [x] Keywords and 2020 Mathematics Subject Classification supplied.
- [x] Software availability, reproducibility, conflict-of-interest, and data-availability statements supplied.
- [x] Figure description supplied for accessibility.
- [x] TAPS package usage kept to ACM-supported packages; packages already loaded by `acmart` are not reloaded explicitly.

## Software component

- [x] Installable Python package in `src/` with `pyproject.toml` metadata.
- [x] Version updated to `0.2.0` in `pyproject.toml`, `CITATION.cff`, and `quadratic_sos.__version__`.
- [x] MIT licence and citation metadata supplied.
- [x] No runtime dependency outside the Python standard library for the core routines.
- [x] Non-interactive smoke-test driver: `python scripts/run_all_checks.py`.
- [x] Smoke-test output separated from deterministic paper data.
- [x] Regression tests and independent brute-force validation scripts supplied.
- [x] Data and figure reproduction scripts supplied.
- [x] Executed notebook supplied for examples.

## Suggested TOMS handling note

The submission is best classified as an Algorithm/software paper: it gives an exact algorithm,
proves correctness, supplies an installable implementation, and includes a reproducibility package
with tests and deterministic computational artefacts.
