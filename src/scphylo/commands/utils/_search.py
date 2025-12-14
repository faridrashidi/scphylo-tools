import os

import click
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm

import scphylo as scp

# from scphylo.pl._trees import _newick_info2_mutation_list


def run_scistree(df_in, alpha, beta, outfile):
    type(outfile)
    # scp.settings.logfile = f"{outfile}/fn_{beta}-fp_{alpha}.log"
    scp.settings.verbosity = "error"
    df_out = scp.tl.scistree(df_in, alpha, beta)
    nll = scp.ul.calc_nll_matrix(df_in, df_out, alpha, beta)
    # scp.io.write(df_out, f"{outfile}/fn_{beta}-fp_{alpha}.CFMatrix")

    # tree = scp.ul.to_tree(df_out)
    # newick, info2, _ = _newick_info2_mutation_list(tree)
    # with open(f"{outfile}/fn_{beta}-fp_{alpha}.newick", "w") as fout:
    #     fout.write(newick + "\n")
    # info2.to_csv(f"{outfile}/fn_{beta}-fp_{alpha}.info2", index=None)
    return nll, df_out, alpha, beta


@click.command(short_help="Grid search for all parameters.")
@click.argument(
    "genotype_file",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
@click.option(
    "--n_threads",
    "-p",
    default=-1,
    type=int,
    show_default=True,
    help="Number of threads.",
)
def search(genotype_file, n_threads):
    """Grid search for all parameters of alpha and beta.

    scphylo search input.SC
    """
    outfile = os.path.splitext(genotype_file)[0]
    # scp.ul.mkdir(outfile)

    df_in = scp.io.read(genotype_file)

    betas = [0.1, 0.2, 0.3, 0.4]
    alphas = [0.1, 0.01, 0.001, 0.0001, 0.00001]
    n_samples = len(betas) * len(alphas)
    if n_threads == -1:
        n_threads = n_samples

    with scp.ul.tqdm_joblib(
        tqdm(
            ascii=True,
            ncols=100,
            desc="GRID SEARCHING:",
            total=n_samples,
            position=0,
        )
    ):
        output = Parallel(n_jobs=n_threads)(
            delayed(run_scistree)(df_in, alpha, beta, outfile)
            for alpha in alphas
            for beta in betas
        )
        output_i = np.argmin([x[0] for x in output])
        df_out = output[output_i][1]
        alpha = output[output_i][2]
        beta = output[output_i][3]
        scp.io.write(df_out, f"{outfile}.scistree.CFMatrix")
        scp.settings.verbosity = "info"
        scp.settings.logfile = f"{outfile}.scistree.log"
        scp.logg.info(f"running ScisTree with alpha={alpha}, beta={beta}, n_threads=1")
        scp.ul.stat(df_in, df_out, alpha, beta, 0)

    return None
