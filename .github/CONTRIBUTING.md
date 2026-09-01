# Contributing guidelines

Welcome to the [scphylo-tools project](https://github.com/faridrashidi/scphylo-tools)!
Before sending your pull requests, please make sure that you **read the full
guidelines**. If you have any questions regarding the contributing guide, please
feel free to [open an issue](https://github.com/faridrashidi/scphylo-tools/issues/new/choose).

## Table of contents

- [Contributing to scphylo-tools](#contributing-to-scphylo-tools)
- [Codebase structure](#codebase-structure)
- [Code style guide](#code-style-guide)
- [Testing](#testing)
- [Writing documentation](#writing-documentation)
- [Building on a new machine](#building-on-a-new-machine)
- [Building C++ files](#building-c-files)
- [Submitting a PR](#submitting-a-pr)
- [Creating a new release](#creating-a-new-release)

## Contributing to scphylo-tools

Clone the scphylo-tools repository from source:

```bash
git clone https://github.com/faridrashidi/scphylo-tools
cd scphylo-tools
```

Install in development mode:

```bash
pixi install --locked
pixi run pre-commit-install
```

## Codebase structure

The scphylo-tools package structure is organized as follows:

- [`src/scphylo`](../src/scphylo): the root of the package.
  - [`src/scphylo/io`](../src/scphylo/io): the input/output module, which
    provides functions for reading and writing data.
  - [`src/scphylo/pl`](../src/scphylo/pl): the plotting module, which provides
    functionality for plotting trees in clonal or dendrogram formats.
  - [`src/scphylo/pp`](../src/scphylo/pp): the preprocessing module, which
    provides functions for filtering and preprocessing data.
  - [`src/scphylo/tl`](../src/scphylo/tl): the tools module, which provides a
    high-level API to compute conflict-free solutions and calculate the
    probability of mutations seeding particular cells.
  - [`src/scphylo/ul`](../src/scphylo/ul): the utilities module, which provides
    utility functions.
  - [`src/scphylo/commands`](../src/scphylo/commands): the CLI commands module,
    which enables running scphylo through the command-line interface.
  - [`src/scphylo/datasets`](../src/scphylo/datasets): the datasets module,
    which contains published single-cell datasets and simulation generators.

Tests structure:

- [`tests`](../tests): the root of the test files.

## Code style guide

We use `ruff` for linting, import sorting, and formatting, and `ty` for static
type checking. Both are integrated as pre-commit hooks. You can use `pre-commit`
to check your changes:

```bash
pixi run lint
```

Please remember to use tags like `TODO:` and `FIXME:` to highlight areas
requiring future attention.

## Testing

We use `pytest` for automated testing. To execute the tests, run:

```bash
pixi run test
```

## Writing documentation

We use NumPy-style docstrings for documentation, with the following
modifications:

- When referring to an argument within the same docstring, enclose that
  reference in double backticks, such as \`\`argument\`\`.
- Prefer adding references to
  [`docs/source/references.bib`](../docs/source/references.bib) instead of
  listing them under the `References` section of the docstring.

To build the documentation, run:

```bash
pixi run docs
```

## Building on a new machine

Pixi installs the compilers and all package dependencies and builds the editable
package from the checked-in C++ sources:

```bash
pixi install --locked
```

To create and validate the wheel and source distribution that are uploaded to
PyPI, run:

```bash
pixi run package
```

## Building C++ files

Normal installs use the checked-in generated files. After changing a `.pyx`
source, regenerate its `.cpp` file explicitly:

```bash
pixi run cythonize
```

## Submitting a PR

Before submitting a new pull request, please make sure you followed these
instructions:

- Make sure that your code follows the conventions above (see
  [Code style guide](#code-style-guide) and
  [Writing documentation](#writing-documentation)).
- If applicable, make sure you have added or modified at least one test to
  account for your changes.
- Make sure that all tests pass locally (see [Testing](#testing)).
- If there is no existing issue that this PR solves, create a new
  [issue](https://github.com/faridrashidi/scphylo-tools/issues/new) briefly
  explaining the problem.

## Creating a new release

Depending on the release type (major, minor, or patch), run:

```bash
pixi run release patch
```

Replace `patch` with `minor` or `major` for larger releases.

By default, this creates a new tag, updates `__version__` wherever necessary,
and commits the changes. After bumping the version, push the commit **and** the
newly created tag upstream. You can do this by setting `push.followTags=true` in
your Git configuration or by using `git push --atomic <branch> <tag>`.
