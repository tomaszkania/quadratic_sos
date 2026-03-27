# quadratic_sos

**Exact integral lengths of sums of squares in real quadratic rings**

This repository contains the reference implementation accompanying the paper
*Computing Integral Length of Sums of Squares and Pythagoras Elements in Real Quadratic Rings*.
It provides exact arithmetic in maximal real quadratic orders, element-wise exact-length
computation, recovery of explicit minimal decompositions, and routines for the integral
Pythagoras number and Pythagoras elements. The implementation uses integer-only row-wise
search in the Minkowski embedding together with unit-square balancing by a positive norm-one
Pell unit. The supplementary table script imports the package directly, so the published
computations and the library use the same exact-arithmetic engine.

## Repository layout

| Path | Purpose |
|---|---|
| `src/quadratic_sos/` | Core package |
| `tests/test_all.py` | Regression tests reproducing the paper's examples and Table 1 |
| `notebooks/paper_illustrations.ipynb` | Executed notebook illustrating the algorithms |
| `scripts/real_quadratic_sos_table.py` | Standalone script reproducing the exact-length table for bounded trace |
| `scripts/balancing_ablation.py` | Pell-orbit balancing ablation from the paper |
| `scripts/cross_field_distributions.py` | Cross-field exact-length distributions for bounded trace |
| `scripts/validation_sweep.py` | Exhaustive bounded-trace validation against brute force |
| `paper/` | LaTeX source and PDF of the manuscript |

## What it computes

For a squarefree integer `D >= 2` and a nonzero element `alpha` of the ring of
integers `O_D` of `Q(sqrt(D))`, the package computes:

- **Representability**: whether `alpha` is a sum of integral squares in `O_D`
- **Exact length**: the minimal `n` such that `alpha = x_1^2 + ... + x_n^2` with `x_i in O_D`
- **Decomposition**: an explicit minimal representation
- **Exact search**: row-wise admissible-parallelogram enumeration without floating-point arithmetic
- **Balancing**: automatic reduction of skewed inputs by even powers of a norm-one Pell unit
- **Pythagoras number**: `P(O_D) in {3,4,5}` and a verified witness element
- **Bulk statistics**: exact-length distributions up to a trace bound

The constructor `RealQuadraticOrder(D)` rejects non-squarefree `D`, and the
functions `exact_length()` and `decomposition()` follow the paper in being
restricted to nonzero inputs.

## Quick start

```python
from quadratic_sos import RealQuadraticOrder, exact_length, decomposition

O = RealQuadraticOrder(10)
alpha = O.elem(18, 2)  # 18 + 2*sqrt(10)

print(exact_length(O, alpha.pair))
# 5

print(decomposition(O, O.from_rational(7).pair))
# e.g. [(2, 0), (1, 0), (1, 0), (1, 0)]
```

## Installation

```bash
pip install -e .
```

## Running tests

```bash
python tests/test_all.py
```

## Reproducing the computational section

Local table over $\mathcal{O}_{10}$:

```bash
python scripts/real_quadratic_sos_table.py --D 10 --trace-bound 40
```

Balancing ablation on the Pell orbit of `18 + 2*sqrt(10)`:

```bash
python scripts/balancing_ablation.py --D 10 --a 18 --b 2 --max-k 3
```

Cross-field distributions:

```bash
python scripts/cross_field_distributions.py --fields 5 10 13 29 53 101 --trace-bound 80
```

Validation sweep against brute force:

```bash
python scripts/validation_sweep.py --max-D 41 --trace-bound 24
```

## Citation

If you use this repository, please cite the accompanying manuscript in `paper/`.
GitHub citation metadata is also provided in `CITATION.cff`. The repository homepage is
<http://github.com/tomaszkania/quadratic_sos>.

## Author

Tomasz Kania  
Institute of Mathematics, Czech Academy of Sciences  
Institute of Mathematics, Jagiellonian University

## License

MIT
