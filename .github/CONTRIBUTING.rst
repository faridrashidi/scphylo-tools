Contributing guidelines
~~~~~~~~~~~~~~~~~~~~~~~

Welcome to the `scphylo-tools project <https://github.com/faridrashidi/scphylo-tools>`_!
Before sending your pull requests, please make sure that you **read the full
guidelines**. If you have any questions regarding the contributing guide, please
feel free to `open an
issue <https://github.com/faridrashidi/scphylo-tools/issues/new/choose>`_.

Table of Contents
=================
- `Contributing to scphylo-tools`_
- `Codebase structure`_
- `Code style guide`_
- `Testing`_
- `Writing documentation`_
- `Building on a new machine`_
- `Building cpp files`_
- `Submitting a PR`_
- `Creating a new release`_


Contributing to scphylo-tools
-----------------------------
Clone the scphylo-tools repository from source::

    git clone https://github.com/faridrashidi/scphylo-tools
    cd scphylo-tools

Install in development mode::

    pixi install --locked
    pixi run pre-commit-install


Codebase structure
------------------
The scphylo-tools package structure is organized as follows:

- `scphylo <scphylo>`_: the root of the package.

    - `scphylo/io <scphylo/io>`_: the input/output module, provides functions
        for reading and writing data.
    - `scphylo/pl <scphylo/pl>`_: the plotting module, provides functionality
        for plotting trees in clonal or dendrogram formats.
    - `scphylo/pp <scphylo/pp>`_: the preprocessing module, provides functions
        for filtering and preprocessing data.
    - `scphylo/tl <scphylo/tl>`_: the tools module, provides a high-level API
        to compute conflict-free solutions and calculate the probability of
        mutations seeding particular cells.
    - `scphylo/ul <scphylo/ul>`_: the utils module, provides utility functions.
    - `scphylo/commands <scphylo/commands>`_: the CLI commands module,
        enables running scphylo via the command-line interface (CLI).
    - `scphylo/datasets <scphylo/datasets>`_: the datasets module, contains
        published single-cell datasets and simulation generators.

Tests structure:

- `tests <tests>`_: the root of the test files.


Code style guide
----------------
We use ``ruff`` for linting, import sorting, and formatting, and ``ty`` for
static type checking. Both are integrated as pre-commit hooks. You can use
``pre-commit`` to check your changes::

    pixi run lint

Please remember to use tags like ``TODO[colon]`` and ``FIXME[colon]`` to
highlight areas requiring future attention.


Testing
-------
We use ``pytest`` for automated testing. To execute the tests, run::

    pixi run test


Writing documentation
---------------------
We use ``numpy``-style docstrings for documentation, with the following
modifications:

- when referring to some argument within the same docstring, enclose that
    reference in double backticks (``).

- prefer adding references to ``references.bib`` instead of listing them under
    the ``References`` section of the docstring.

In order to build the documentation, run::

    pixi run docs


Building on a new machine
-------------------------
Pixi installs the compilers and all package dependencies and builds the editable
package from the checked-in C++ sources::

    pixi install --locked

To create and validate the wheel and source distribution that are uploaded to
PyPI, run::

    pixi run package


Building cpp files
------------------
Normal installs use the checked-in generated files. After changing a ``.pyx``
source, regenerate its ``.cpp`` file explicitly::

    pixi run cythonize


Submitting a PR
---------------
Before submitting a new pull request, please make sure you followed these
instructions:

- make sure that your code follows the above specified conventions (see
    `Code style guide`_ and `Writing documentation`_).
- if applicable, make sure you've added/modified at least 1 test to account
    for the changes you've made.
- make sure that all tests pass locally (see `Testing`_).
- if there is no existing issue that this PR solves, create a new
    `issue <https://github.com/faridrashidi/scphylo-tools/issues/new>`_ briefly
    explaining the problem.


Creating a new release
----------------------
Depending on the release type (major, minor, or patch), run::

    pixi run release patch

Replace ``patch`` with ``minor`` or ``major`` for larger releases.

By default, this will create a new tag and automatically update the
``__version__`` wherever necessary, commit the changes and create a new tag.
After bumping the version, push the commit **AND** the newly created tag to
upstream. This can be done by setting ``push.followtags=true`` in your git
config or by using ``git push --atomic <branch> <tag>``.
