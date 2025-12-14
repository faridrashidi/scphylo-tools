CLI
===

The scPhylo-tools package includes a command-line interface (CLI).
Once scphylo is installed (see :ref:`installation tutorial <installationguide>`),
the ``scphylo`` command becomes available in your terminal. ``scphylo`` is a
command-line tool featuring several sub-commands. You can view a list of all
available commands by typing ``scphylo --help``. This produces the
following output:


.. click:: scphylo.commands.scphylo:cli
    :prog: scphylo
    :commands: --help


``mcalling`` - Run Mutation Calling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. click:: scphylo.commands.scphylo:cli
    :prog: scphylo
    :commands: mcalling
    :nested: full
