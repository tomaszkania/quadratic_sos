.PHONY: install test reproduce validate smoke paper clean

PAPER = integral_sos_real_quadratic

install:
	python -m pip install -e .[repro]

test:
	python -m pytest -q

reproduce:
	python scripts/reproduce_tables.py

validate:
	python scripts/validation_sweep.py --max-D 41 --trace-bound 24

smoke:
	python scripts/run_all_checks.py

paper:
	cd paper && pdflatex -interaction=nonstopmode $(PAPER).tex && (bibtex $(PAPER) || bibtex8 $(PAPER)) && pdflatex -interaction=nonstopmode $(PAPER).tex && pdflatex -interaction=nonstopmode $(PAPER).tex

clean:
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +
	find . -name '*.pyc' -delete
	rm -rf .pytest_cache src/*.egg-info build dist data/smoke
	rm -f paper/*.aux paper/*.log paper/*.out paper/*.toc paper/*.bbl paper/*.blg paper/*.fdb_latexmk paper/*.fls
