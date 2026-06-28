VALID_RELEASE_PARTS := patch minor major
RELEASE_ARGS = $(filter-out $@,$(MAKECMDGOALS))

.PHONY: help clean lint test doc build install install-dev cythonize release $(VALID_RELEASE_PARTS)

help:
	@echo "available commands"
	@echo " - lint         : run linting and flaking"
	@echo " - test         : run all unit tests"
	@echo " - doc          : build the documentation"
	@echo " - release      : bump version and create a release tag"
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
	rm -rf docs/source/sg_execution_times.rst

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
	rm -rf docs/source/sg_execution_times.rst
	cd docs && $(MAKE) clean html
	cd docs/build/html && python -m http.server 8080

build:
	python setup.py sdist

install:
	python -m pip install .

install-dev:
	git fetch
	python -m pip install -e '.[dev,docs]'
	pre-commit install

cythonize:
	CYTHONIZE=1 pip install -e .

release:
	@if [ "$(words $(RELEASE_ARGS))" -ne 1 ] || [ -n "$(filter-out $(VALID_RELEASE_PARTS),$(RELEASE_ARGS))" ]; then \
		echo "usage: make release [patch|minor|major]"; \
		exit 1; \
	fi
	uv run --extra dev bump-my-version bump $(RELEASE_ARGS)

$(VALID_RELEASE_PARTS):
	@:
