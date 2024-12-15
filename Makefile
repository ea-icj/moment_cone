.PHONY: all_tests unittest doctest mypy fixme todo doc

all_tests: unittest doctest mypy fixme todo


unittest:
	python -m unittest discover --verbose .

doctest:
	python -m pytest --doctest-modules --verbose src/cone

mypy:
	python -m mypy src/cone tests

fixme:
	grep -r --exclude-dir=__pycache__ --color --line-number "FIXME" tests/ src/

todo:
	grep -r --exclude-dir=__pycache__ --color --line-number "TODO" tests/ src/

doc:
	pdoc3 --html --force src/cone
