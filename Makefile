PDFLATEX ?= pdflatex
BIBTEX ?= bibtex
PAPER = exact_diagonal_representability_real_quadratic

.PHONY: install test reproduce validate smoke paper submission clean

install:
	python -m pip install -e .[repro]

test:
	python -m pytest -q

reproduce:
	python scripts/reproduce_tables.py

validate:
	python scripts/validation_sweep.py

smoke:
	python scripts/run_all_checks.py

paper:
	cd paper && $(PDFLATEX) -interaction=nonstopmode -halt-on-error $(PAPER).tex
	cd paper && ($(BIBTEX) $(PAPER) || bibtex.original $(PAPER))
	cd paper && $(PDFLATEX) -interaction=nonstopmode -halt-on-error $(PAPER).tex
	cd paper && $(PDFLATEX) -interaction=nonstopmode -halt-on-error $(PAPER).tex

submission:
	python scripts/make_submission_archive.py

clean:
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +
	find . -name '*.pyc' -delete
	find . -name '*.pyo' -delete
	rm -rf .pytest_cache build dist src/*.egg-info
	rm -f paper/*.aux paper/*.bbl paper/*.log paper/*.out paper/*.toc paper/*.blg paper/*.fls paper/*.fdb_latexmk
