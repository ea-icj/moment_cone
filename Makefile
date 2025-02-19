.PHONY: all_tests unittest doctest mypy fixme todo doc

all_tests: unittest doctest mypy fixme todo


unittest:
	python3 -m unittest discover --verbose .

doctest:
	python3 -m pytest --doctest-modules --verbose --ignore=src/cone/main_to_be_inserted.py src/cone

mypy:
	python3 -m mypy src/cone tests

fixme:
	grep -ir --exclude-dir=__pycache__ --color --line-number "FIXME" tests/ src/

todo:
	grep -ir --exclude-dir=__pycache__ --color --line-number "TODO" tests/ src/

doc:
	pdoc3 --html --force src/cone

tout_test:
	python -m pytest --doctest-modules --verbose --ignore=src/tout_test/reference_datas src/tout_test
	python -m mypy --check-untyped-defs src/tout_test
