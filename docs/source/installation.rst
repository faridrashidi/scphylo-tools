:tocdepth: 1

Installation
------------

scPhylo-tools requires Python 3.7 or later.

Using PyPI
^^^^^^^^^^

Install scphylo from PyPI_ using::

    pip install -U scphylo

``-U`` is short for ``--upgrade``.
If you get a ``Permission denied`` error, use
``pip install -U scphylo --user`` instead.


Development Version
^^^^^^^^^^^^^^^^^^^

To work with the latest development version, install from GitHub_ using::

    git clone https://github.com/faridrashidi/scphylo
    cd scphylo
    make install-dev

``-e`` stands for ``--editable`` and makes sure that your environment
is updated when you pull new changes from GitHub. The ``'[dev]'`` options
installs requirements needed for development, because scPhylo-tools
is bundled with an additional library.

Your contributions to improve scPhylo-tools is highly appreciated! Please
check out our `contributing guide`_.

If you run into issues, do not hesitate to approach us or
raise a `GitHub issue`_.

.. _PyPI: https://pypi.org/project/scphylo
.. _Github: https://github.com/faridrashidi/scphylo
.. _Github issue: https://github.com/faridrashidi/scphylo/issues/new/choose
.. _contributing guide: https://github.com/faridrashidi/scphylo/blob/master/CONTRIBUTING.rst
