.. _installationguide:

Installation
============

scPhylo-tools requires Python 3.11 or later.

Using PyPI
^^^^^^^^^^

To install scphylo-tools from PyPI_, run::

    pip install scphylo-tools

If you encounter a ``Permission denied`` error, use
``pip install scphylo-tools --user`` instead.

The PyPI package is suitable for the Python-only functionality. Native Graphviz
and the R packages used by ``clonal_tree`` and ``dendro_tree`` cannot be
installed by pip alone; use the locked Pixi environment below when you need the
complete plotting stack.


Using Pixi
^^^^^^^^^^

Pixi_ is the recommended way to develop scPhylo-tools and to use all plotting
features. After installing Pixi, clone the repository and install its lockfile::

    git clone https://github.com/faridrashidi/scphylo-tools
    cd scphylo-tools
    pixi install --locked

The default environment includes Python, Graphviz, ``pygraphviz``, R,
``rpy2``, ``ggtree``, ``ggtreeExtra``, ``graph-tool``, and MPI. Verify the
complete development setup with::

    pixi run test

The locked Pixi environments support Linux and macOS. Windows users can use WSL
for the complete native and R plotting stack.

Development Version
^^^^^^^^^^^^^^^^^^^

To work with the latest development version, use the Pixi setup above and then
install the pre-commit hooks::

    pixi run pre-commit-install

Your contributions to improve scPhylo-tools are highly appreciated! Please
refer to our `contributing guide`_.

If you encounter any issues, please feel free to open a `GitHub issue`_.

.. _PyPI: https://pypi.org/project/scphylo-tools
.. _Pixi: https://pixi.prefix.dev/latest/installation/
.. _Github: https://github.com/faridrashidi/scphylo-tools
.. _Github issue: https://github.com/faridrashidi/scphylo-tools/issues/new/choose
.. _contributing guide: https://github.com/faridrashidi/scphylo-tools/blob/main/.github/CONTRIBUTING.rst
