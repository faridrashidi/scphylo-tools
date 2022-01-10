from collections import OrderedDict

import click

import scphylo as scp
from scphylo.commands.caller import cli_caller
from scphylo.commands.methyl import cli_methyl
from scphylo.commands.solver import cli_solver
from scphylo.commands.utils import cli_utils


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


@click.version_option(version=scp.__version__)
@click.group(
    cls=NaturalOrderGroup,
    commands=OrderedDict(),
    context_settings={"max_content_width": 300, "terminal_width": 300},
)
def cli():
    """scPhylo-tools.

    A python toolkit for single-cell tumor phylogenetic analysis
    """
    return None


cli.add_command(cli_solver)
cli.add_command(cli_caller)
cli.add_command(cli_utils)
cli.add_command(cli_methyl)
