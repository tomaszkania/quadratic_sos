# Commit name

toms-readiness-v0.2.0

# Description

Convert the manuscript and companion implementation into a TOMS-ready submission bundle. The manuscript now uses the ACM `acmart` primary article class with TOMS metadata, CCS descriptors, keywords, accessibility description for the figure, ACM reference formatting, and explicit software/reproducibility statements. The Python package version is updated to `0.2.0`, includes typed package metadata, and has refreshed documentation, artefact instructions, a submission checklist, and citation metadata.

The implementation update fixes the validation workflow so the brute-force checker is independent of Peters' representability criterion, the production row-search engine, and the closed-form Pythagoras-number classification. Smoke-test outputs are separated from deterministic paper data, quick smoke tests no longer overwrite the paper figure, and the test suite and reproduction scripts are updated accordingly.
