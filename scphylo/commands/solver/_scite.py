import os

import click
import numpy as np
from joblib import Parallel, delayed

import scphylo as scp


@click.command(short_help="Run SCITE.")
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
    "--n_iters",
    "-l",
    default=1000000,
    type=int,
    show_default=True,
    help="Number of iterations.",
)
@click.option(
    "--n_restarts",
    "-r",
    default=3,
    type=int,
    show_default=True,
    help="Number of restarts.",
)
@click.option(
    "--experiment",
    "-e",
    is_flag=True,
    default=False,
    type=bool,
    show_default=True,
    help="Is in experiment mode.",
)
@click.option(
    "--time_limit",
    "-t",
    default=86400,
    type=float,
    show_default=True,
    help="Time limit for the experiment part (in seconds).",
)
@click.option(
    "--smooth_rate",
    "-s",
    default=2,
    type=float,
    show_default=True,
    help="Smooth rate for the experiment part.",
)
@click.option(
    "--iters_rate",
    "-i",
    default=30000,
    type=int,
    show_default=True,
    help="Interations rate for the experiment part.",
)
def scite(
    genotype_file,
    alpha,
    beta,
    n_iters,
    n_restarts,
    experiment,
    time_limit,
    smooth_rate,
    iters_rate,
):
    """SCITE.

    Tree inference for single-cell data :cite:`SCITE`.

    scphylo scite input.SC 0.0001 0.1 -l 1000000 -r 3 -e -t 86400 -s 2
    """
    outfile = os.path.splitext(genotype_file)[0]

    scp.settings.verbosity = "info"

    df_in = scp.io.read(genotype_file)
    if not experiment:
        scp.settings.logfile = f"{outfile}.scite.log"
        df_out = scp.tl.scite(
            df_in,
            alpha=alpha,
            beta=beta,
            n_iters=n_iters,
            n_restarts=n_restarts,
        )
        scp.io.write(df_out, f"{outfile}.scite.CFMatrix")
    else:
        scp.settings.logfile = f"{outfile}.scite.log"
        df_out, running_time, _, _ = scp.tl.scite(
            df_in,
            alpha=alpha,
            beta=beta,
            n_iters=iters_rate,
            n_restarts=1,
            experiment=True,
        )
        n_iters = int(smooth_rate * iters_rate * time_limit / running_time)

        def run(i):
            do, r, s, b = scp.tl.scite(
                df_in,
                alpha=alpha,
                beta=beta,
                n_iters=n_iters,
                n_restarts=1,
                experiment=True,
            )
            return do, r, s, b

        output = Parallel(n_jobs=n_restarts)(delayed(run)(i) for i in range(n_restarts))

        scores = [x[2] for x in output]
        betas = [x[3] for x in output]
        best_i = np.argmax(scores)
        df_out = output[best_i][0]

        scp.ul.stat(df_in, df_out, alpha, beta, output[best_i][1])
        scp.logg.info(f"score: {output[best_i][2]}")
        scp.logg.info(f"beta: {output[best_i][3]}")
        scp.logg.info(f"n_iters: {n_iters}")
        scp.logg.info(f"scores: {','.join(list(map(str, scores)))}")
        scp.logg.info(f"betas: {','.join(list(map(str, betas)))}")
        scp.logg.info(f"picked: {best_i}")

        scp.io.write(df_out, f"{outfile}.scite.CFMatrix")

    return None
