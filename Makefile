help:
	@echo "available commands"
	@echo " - lint         : run linting and flaking"
	@echo " - test         : run all unit tests"
	@echo " - doc          : build the documentation"
	@echo " - install      : install the package"
	@echo " - install-dev  : install the package in dev mode"
	@echo " - cythonize    : cythonize all pyx files"
	@echo " - clean        : clean the repo"

clean:
	rm -rf docs/build
	rm -rf htmlcov
	rm -rf .coverage*
	rm -rf docs/source/scphylo*
	rm -rf docs/source/auto_examples
	rm -rf docs/source/gen_modules

lint:
	pre-commit run --all-files

test:
	rm -rf htmlcov
	pytest --disable-warnings --cov=scphylo --cov-report=html ./tests
	rm -rf .coverage*

doc:
	rm -rf docs/source/scphylo*
	rm -rf docs/source/auto_examples
	rm -rf docs/source/gen_modules
	cd docs && $(MAKE) clean html
	cd docs/build/html && python -m http.server 8080

install:
	python -m pip install .

install-dev:
	git fetch
	python -m pip install -e '.[dev,docs]'
	pre-commit install

cythonize:
	CYTHONIZE=1 pip install -e .
