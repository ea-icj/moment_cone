.PHONY: all_tests unittest doctest mypy fixme todo doc

all_tests: unittest doctest mypy fixme todo


unittest:
	python3 -m unittest discover --verbose .

doctest:
	python3 -m pytest --doctest-modules --verbose --ignore=src/cone/reference_datas src/cone

mypy:
	python3 -m mypy --check-untyped-defs --exclude=src/cone/reference_datas src/cone tests

mypy_strict:
	python3 -m mypy --strict --no-warn-unused-ignores src/cone tests


fixme:
	grep -ir --exclude-dir=__pycache__ --color --line-number "FIXME" tests/ src/

todo:
	grep -ir --exclude-dir=__pycache__ --color --line-number "TODO" tests/ src/

doc:
	pdoc3 --html --force src/cone

