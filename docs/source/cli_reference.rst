CLI
===

A command line interface (CLI) is available in scPhylo-tools package.
After you have scphylo correctly installed on your machine
(see :ref:`installation tutorial <installationguide>`), the ``scphylo``
command will become available in the terminal. ``scphylo`` is a
command line tool with some sub-commands. You can get quick info on all the
available commands typing ``scphylo --help``. You will get the
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
