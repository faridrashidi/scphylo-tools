import os

import click

import scphylo as scp


@click.command(short_help="Run HUNTRESS.")
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
@click.option(
    "--n_threads",
    "-p",
    default=1,
    type=int,
    show_default=True,
    help="Number of threads.",
)
def huntress(genotype_file, alpha, beta, n_threads):
    """HUNTRESS.

    scphylo huntress input.SC 0.0001 0.1 -p 8
    """
    outfile = os.path.splitext(genotype_file)[0]

    scp.settings.verbosity = "info"
    scp.settings.logfile = f"{outfile}.huntress.log"

    df_in = scp.io.read(genotype_file)
    df_out = scp.tl.huntress(df_in, alpha=alpha, beta=beta, n_threads=n_threads)
    scp.io.write(df_out, f"{outfile}.huntress.CFMatrix")

    return None
