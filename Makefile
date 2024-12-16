all_tests: unittest doctest mypy


unittest:
	python -m unittest discover --verbose .

doctest:
	python -m pytest --doctest-modules --verbose --ignore=src/cone/main_to_be_inserted.py src/cone

mypy:
	python -m mypy src/cone tests

doc:
	pdoc3 --html --force src/cone
