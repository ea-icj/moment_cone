.PHONY: all_tests unittest doctest mypy fixme todo doc

all_tests: unittest doctest mypy fixme todo


unittest:
	python3 -m unittest discover --verbose .

doctest:
	python3 -m pytest --doctest-modules --verbose --ignore=src/moment_cone/reference_datas src/moment_cone

mypy:
	python3 -m mypy --install-types --check-untyped-defs --exclude=src/moment_cone/reference_datas/ineq_ src/moment_cone tests

mypy_strict:
	python3 -m mypy --install-types --strict --no-warn-unused-ignores --exclude src/moment_cone/reference_datas/ineq_ src/moment_cone tests


fixme:
	grep -ir --exclude-dir=__pycache__ --color --line-number "FIXME" tests/ src/

todo:
	grep -ir --exclude-dir=__pycache__ --color --line-number "TODO" tests/ src/

doc:
	pdoc3 --html --force src/moment_cone

