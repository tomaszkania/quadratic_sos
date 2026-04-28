.PHONY: install test reproduce validate smoke paper clean

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
	cd paper && pdflatex -interaction=nonstopmode *.tex && pdflatex -interaction=nonstopmode *.tex

clean:
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +
	find . -name '*.pyc' -delete
	rm -f paper/*.aux paper/*.log paper/*.out paper/*.toc
