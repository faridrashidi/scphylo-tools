import os

import click

import scphylo as scp


@click.command(short_help="Run OncoNEM.")
@click.argument(
    "genotype_file",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument(
    "alpha",
    required=True,
    type=float,
)
@click.argument(
    "beta",
    required=True,
    type=float,
)
def onconem(genotype_file, alpha, beta):
    """OncoNEM.

    Inferring tumor evolution from single-cell sequencing data :cite:`OncoNEM`.

    scphylo onconem input.SC 0.0001 0.1
    """
    outfile = os.path.splitext(genotype_file)[0]

    scp.settings.verbosity = "info"
    scp.settings.logfile = f"{outfile}.onconem.log"

    df_in = scp.io.read(genotype_file)
    df_out = scp.tl.onconem(df_in, alpha=alpha, beta=beta)
    scp.io.write(df_out, f"{outfile}.onconem.CFMatrix")

    return None
