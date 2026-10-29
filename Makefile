.PHONY: install test reproduce verify
install:
	pip install -r requirements.txt
test:
	python -m unittest discover -s tests -v
reproduce:
	python scripts/run_analysis.py
verify:
	python scripts/verify.py
