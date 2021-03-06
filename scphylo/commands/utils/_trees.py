import os

import click

import scphylo as scp
from scphylo.pl._trees import _newick_info2_mutation_list


@click.command(short_help="Convert conflict-free to newick file.")
@click.argument(
    "cf_file",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
def cf2newick(cf_file):
    """Convert conflict-free matrix to newick."""
    outfile = os.path.splitext(cf_file)[0]
    cf = scp.io.read(cf_file)
    tree = scp.ul.to_tree(cf)
    newick, info2, mutations = _newick_info2_mutation_list(tree)
    with open(f"{outfile}.newick", "w") as fout:
        fout.write(newick + "\n")
    info2.to_csv(f"{outfile}.info2", index=None)

    return None


@click.command(short_help="Convert conflict-free to clonal tree.")
@click.argument(
    "cf_file",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
def cf2tree(cf_file):
    """Draw conflict-free matrix to clonal tree."""
    outfile = os.path.splitext(cf_file)[0]
    cf = scp.io.read(cf_file)
    tree = scp.ul.to_tree(cf)
    scp.pl.clonal_tree(tree, output_file=f"{outfile}.png")

    return None
