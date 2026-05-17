# quadratic_sos

**Exact integral lengths of sums of squares in real quadratic rings**

`quadratic_sos` is the reference implementation accompanying the TOMS submission
*Computing Integral Length of Sums of Squares and Pythagoras Elements in Real Quadratic Rings*.
It provides exact arithmetic in maximal real quadratic orders, element-wise exact-length
computation, recovery of explicit minimal decompositions, and routines for integral Pythagoras
numbers and Pythagoras elements.

The implementation is deliberately conservative: the core package has no runtime dependency outside
the Python standard library, all comparisons are performed with integer arithmetic, and the public
constructor accepts only squarefree radicands `D >= 2`.

## Version

Current repository version: **0.2.0**.

This version is the TOMS-readiness update. It aligns the manuscript with the ACM primary article
template, separates deterministic paper data from smoke-test output, and strengthens the independent
brute-force validation path.

## Repository layout

| Path | Purpose |
|---|---|
| `src/quadratic_sos/` | Installable core package |
| `tests/test_all.py` | Regression tests for examples, witnesses, and table counts |
| `scripts/run_all_checks.py` | Non-interactive TOMS smoke-test driver |
| `scripts/reproduce_tables.py` | Regenerates deterministic table data and the balancing plot |
| `scripts/validation_sweep.py` | Brute-force bounded-trace validation driver |
| `notebooks/paper_illustrations.ipynb` | Executed notebook illustrating the algorithms |
| `data/` | Deterministic paper data and smoke-test summaries |
| `paper/` | ACM/TOMS LaTeX source, figure files, and compiled PDF |

## What it computes

For a squarefree integer `D >= 2` and a nonzero element `alpha` of the ring of integers `O_D` of
`Q(sqrt(D))`, the package computes:

- representability of `alpha` as a sum of integral squares in `O_D`;
- the exact length, i.e. the minimal `n` with `alpha = x_1^2 + ... + x_n^2`;
- an explicit minimal decomposition when one exists;
- exact row-wise enumeration of the admissible search region in the Minkowski embedding;
- unit-square balancing using a positive norm-one Pell unit;
- `P(O_D) in {3, 4, 5}` and a verified Pythagoras element;
- bounded-trace exact-length distributions used in the computational section of the paper.

The functions `exact_length()` and `decomposition()` follow the paper and reject the zero element.
Non-representability is returned as `None`, corresponding to length infinity in the paper.

## Quick start

```python
from quadratic_sos import RealQuadraticOrder, decomposition, exact_length

order = RealQuadraticOrder(10)
alpha = order.elem(18, 2)  # 18 + 2*sqrt(10)

print(exact_length(order, alpha.pair))
# 5

print(decomposition(order, order.from_rational(7).pair))
# [(2, 0), (1, 0), (1, 0), (1, 0)]
```

## Installation

```bash
python -m pip install -e .
```

For tests, figures, and notebook execution, install the reproducibility extras:

```bash
python -m pip install -e .[repro]
```

## Running tests

```bash
python -m pytest -q
```

A successful run currently reports 26 passing tests.

## TOMS smoke test

The referee-facing non-interactive command is:

```bash
python -m pip install -e .[repro]
python scripts/run_all_checks.py
```

The smoke test runs the regression suite, regenerates quick deterministic data under `data/smoke/`,
performs a small brute-force validation sweep, writes `data/run_all_checks_summary.txt`, and ends
with an `OK` line.

## Reproducing the computational section

Regenerate the deterministic data used in the paper:

```bash
python scripts/reproduce_tables.py
```

This writes:

- `data/sos_D10_distribution.txt`;
- `data/sos_balancing_ablation.tsv`;
- `data/sos_cross_field.tsv`;
- `data/sos_validation_summary.tsv`;
- `paper/figures/sos_balancing_rows.pdf` and `.png` when `matplotlib` is available.

Individual drivers are also available:

```bash
python scripts/real_quadratic_sos_table.py --D 10 --trace-bound 40
python scripts/balancing_ablation.py --D 10 --a 18 --b 2 --max-k 3
python scripts/cross_field_distributions.py --fields 5 10 13 29 53 101 --trace-bound 80
python scripts/validation_sweep.py --max-D 41 --trace-bound 24
```

Optional SageMath and PARI/GP probes are intentionally outside the mandatory smoke-test path:

```bash
python scripts/optional_cas_baselines.py
```

## Building the paper

The manuscript uses the ACM primary article class in manuscript mode for TOMS submission.
From the repository root:

```bash
make paper
```

The compiled PDF is written to `paper/integral_sos_real_quadratic.pdf`.

## Citation

If you use this repository, cite the accompanying manuscript in `paper/`. GitHub citation metadata
is provided in `CITATION.cff`.

## Author

Tomasz Kania  
Institute of Mathematics, Czech Academy of Sciences  
Institute of Mathematics, Jagiellonian University

## Licence

MIT
