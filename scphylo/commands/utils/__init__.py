"""Utils Module."""

from collections import OrderedDict

import click

from scphylo.commands.utils._consensus import consensus
from scphylo.commands.utils._partf import partf
from scphylo.commands.utils._score import score
from scphylo.commands.utils._search import search
from scphylo.commands.utils._trees import cf2newick, cf2tree


class NaturalOrderGroup(click.Group):
    """Command group trying to list subcommands in the order they were added.

    Make sure you initialize the `self.commands` with OrderedDict instance.
    With decorator, use::
        @click.group(cls=NaturalOrderGroup, commands=OrderedDict())
    """

    def list_commands(self, ctx):
        """List command names as they are in commands dict.

        If the dict is OrderedDict, it will preserve the order commands
        were added.
        """
        return self.commands.keys()


@click.group(
    name="utils",
    cls=NaturalOrderGroup,
    commands=OrderedDict(),
    context_settings={"max_content_width": 300, "terminal_width": 300},
)
def cli_utils():
    """Utils module."""
    return None


cli_utils.add_command(score)
cli_utils.add_command(search)
cli_utils.add_command(cf2newick)
cli_utils.add_command(cf2tree)
cli_utils.add_command(consensus)
cli_utils.add_command(partf)
