"""Solver Module."""

from collections import OrderedDict

import click

from scphylo.commands.solver._bnb import bnb
from scphylo.commands.solver._booster import booster
from scphylo.commands.solver._consensus import consensus
from scphylo.commands.solver._gpps import gpps
from scphylo.commands.solver._grmt import grmt
from scphylo.commands.solver._huntress import huntress
from scphylo.commands.solver._onconem import onconem
from scphylo.commands.solver._partf import partf
from scphylo.commands.solver._phiscs import phiscsb, phiscsi
from scphylo.commands.solver._scistree import scistree
from scphylo.commands.solver._scite import scite
from scphylo.commands.solver._score import score
from scphylo.commands.solver._search import search
from scphylo.commands.solver._sphyr import sphyr
from scphylo.commands.solver._trees import cf2newick, cf2tree


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
    name="solver",
    cls=NaturalOrderGroup,
    commands=OrderedDict(),
    context_settings={"max_content_width": 300, "terminal_width": 300},
)
def cli_solver():
    """Solver module."""
    return None


cli_solver.add_command(booster)
cli_solver.add_command(partf)
cli_solver.add_command(phiscsb)
cli_solver.add_command(phiscsi)
cli_solver.add_command(scite)
cli_solver.add_command(scistree)
cli_solver.add_command(bnb)
cli_solver.add_command(huntress)
cli_solver.add_command(grmt)
cli_solver.add_command(sphyr)
cli_solver.add_command(gpps)
cli_solver.add_command(onconem)
cli_solver.add_command(cf2newick)
cli_solver.add_command(cf2tree)
cli_solver.add_command(score)
cli_solver.add_command(consensus)
cli_solver.add_command(search)
