import os

import click

import scphylo as scp


@click.command(short_help="Run PhISCS-BnB.")
@click.argument(
    "genotype_file",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
@click.option(
    "--bounding",
    "-b",
    default="real",
    type=click.Choice(["real", "simulated"]),
    show_default=True,
    help="Number of threads.",
)
@click.option(
    "--time_limit",
    "-t",
    default=None,
    type=int,
    show_default=True,
    help="Time limit.",
)
def bnb(genotype_file, bounding, time_limit):
    """PhISCS-BnB.

    A fast branch and bound algorithm for
    the perfect tumor phylogeny reconstruction problem :cite:`PhISCS-BnB`.

    scphylo bnb input.SC -b real
    """
    outfile = os.path.splitext(genotype_file)[0]

    scp.settings.verbosity = "info"
    scp.settings.logfile = f"{outfile}.bnb.log"

    df_in = scp.io.read(genotype_file)
    df_out = scp.tl.bnb(df_in, bounding=bounding, time_limit=time_limit)
    scp.io.write(df_out, f"{outfile}.bnb.CFMatrix")

    return None
