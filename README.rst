|Tests| |Docs| |PyPI| |Python Version| |License|


scphylo-tools: a python toolkit for single-cell tumor phylogenetic analysis
===========================================================================

.. image:: https://raw.githubusercontent.com/faridrashidi/scphylo-tools/main/docs/source/_static/images/overview.png
   :target: https://github.com/faridrashidi/scphylo-tools
   :alt: overview
   :align: center

**scphylo-tools** is a comprehensive Python toolkit designed for single-cell tumor
phylogenetic analysis. It integrates efficient implementations of state-of-the-art
methods to infer evolutionary trees from single-cell genotype data.

Features
--------

*   **Phylogenetic Inference**: Wrappers and implementations for state-of-the-art
    inference methods.

*   **API & CLI**: Provides a Python API as well as command-line tools.

*   **Visualization**: Built-in visualization tools.

*   **Datasets**: Access to available datasets.

*   **Ease of Use**: A unified Python interface for complex phylogenetic workflows.

Installation
------------

Install the Python package from PyPI:

.. code-block:: bash

    pip install scphylo-tools

For development and the complete native stack, use `Pixi
<https://pixi.prefix.dev/latest/installation/>`_. Pixi installs Python together
with Graphviz, R, ``rpy2``, ``ggtree``, ``ggtreeExtra``, ``graph-tool``, and
MPI from the lockfile:

.. code-block:: bash

    git clone https://github.com/faridrashidi/scphylo-tools
    cd scphylo-tools
    pixi install --locked
    pixi run test

The locked Pixi environments support Linux and macOS. On Windows, use WSL for
the complete native and R plotting stack.

The regular PyPI wheel and source distribution remain available for users who
only need the Python package. See the `installation guide
<https://scphylo-tools.readthedocs.io/en/latest/installation.html>`_ for the
complete development environment and plotting details.

Documentation
-------------

Detailed documentation and tutorials are available at `Read the Docs <https://scphylo-tools.readthedocs.io>`_.

.. |Tests| image:: https://img.shields.io/github/actions/workflow/status/faridrashidi/scphylo-tools/ci.yml?branch=main&logo=github&logoColor=white&style=flat-square&label=tests&labelColor=000000&cacheSeconds=0
    :target: https://github.com/faridrashidi/scphylo-tools/actions/workflows/ci.yml
    :alt: Tests

.. |Docs| image:: https://img.shields.io/website?url=https%3A%2F%2Fscphylo-tools.readthedocs.io%2F&up_message=online&down_message=offline&logo=readthedocs&logoColor=white&style=flat-square&label=docs&labelColor=000000&cacheSeconds=0
    :target: https://scphylo-tools.readthedocs.io/
    :alt: Docs

.. |PyPI| image:: https://img.shields.io/pypi/v/scphylo-tools?logo=pypi&logoColor=white&style=flat-square&labelColor=000000&cacheSeconds=0
    :target: https://pypi.org/project/scphylo-tools/
    :alt: PyPI

.. |Python Version| image:: https://img.shields.io/pypi/pyversions/scphylo-tools?logo=python&logoColor=white&style=flat-square&labelColor=000000&cacheSeconds=0
    :target: https://pypi.org/project/scphylo-tools/
    :alt: Python Version

.. |License| image:: https://img.shields.io/pypi/l/scphylo-tools?logo=creativecommons&logoColor=white&style=flat-square&labelColor=000000&color=blueviolet&cacheSeconds=0
    :target: https://github.com/faridrashidi/scphylo-tools/blob/main/LICENSE.md
    :alt: License
