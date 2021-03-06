"""Huntress Module."""

import pandas as pd

import scphylo as scp
from scphylo.tl.solver.huntress._huntress import Reconstruct


def huntress(df_input, alpha, beta, n_threads=1):
    """Solving using HUNTRESS.

    HUNTRESS: Provably fast intratumor heterogeneity inference from single-cell
    sequencing data

    Parameters
    ----------
    df_input : :class:`pandas.DataFrame`
        Input genotype matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1), absence (0) and missing
        entires (3).
    alpha : :obj:`float`
        False positive error rate.
    beta : :obj:`float`
        False negative error rate.
    n_threads : :obj:`int`, optional
        Number of threads, by default 1

    Returns
    -------
    :class:`pandas.DataFrame`
        A conflict-free matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1) and absence (0).
    """
    scp.logg.info(
        f"running HUNTRESS with alpha={alpha}, beta={beta}, n_threads={n_threads}"
    )

    running_time = 0
    if alpha == 0:
        matrix_out, running_time = Reconstruct(
            df_input,
            Algchoice="FN",
            n_proc=n_threads,
        )
    else:
        fn_conorm = 0.1
        fp_conorm = fn_conorm * alpha / beta
        matrix_out, running_time = Reconstruct(
            df_input,
            n_proc=n_threads,
            fnfp=51,
            post_fn=fn_conorm,
            post_fp=fp_conorm,
        )

    df_output = pd.DataFrame(matrix_out)
    df_output.columns = df_input.columns
    df_output.index = df_input.index
    df_output.index.name = "cellIDxmutID"

    scp.ul.stat(df_input, df_output, alpha, beta, running_time)

    return df_output
