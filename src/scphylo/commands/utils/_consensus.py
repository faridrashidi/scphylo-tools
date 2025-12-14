import click

import scphylo as scp


@click.command(
    short_help="Build consensus tree between two phylogenetic trees (Trisicell-Cons)."
)
@click.argument(
    "first_tree",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument(
    "second_tree",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument(
    "consensus_path",
    required=True,
    type=click.Path(exists=False),
)
def consensus(first_tree, second_tree, consensus_path):
    """Build consensus tree between two phylogenetic trees.

    It writes the conflict-free matrix representing the consensus tree
    into the `consensus_path` filepath.

    scphylo consensus first_tree.CFMatrix second_tree.CFMatrix
    """
    scp.settings.verbosity = "info"

    sc1 = scp.io.read(first_tree)
    sc2 = scp.io.read(second_tree)
    final_tree = scp.tl.consensus(sc1, sc2)
    scp.io.write(final_tree.graph["data"], consensus_path)

    return None
