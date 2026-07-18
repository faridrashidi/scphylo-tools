import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd

import scphylo as scp


def siclonefit(
    df_input,
    alpha,
    beta,
    n_restarts=3,
    n_iters=500,
    n_burnin=100,
    return_tree=False,
    experiment=False,
    jar_path=None,
    java_executable=None,
    tools_dir=None,
):
    """Solving using SiCloneFit.

    Bayesian inference of population structure, genotype, and phylogeny of tumor clones
    from single-cell genome sequencing data :cite:`SiCloneFit`.

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
    n_restarts : :obj:`int`, optional
        Number of restarts, by default 3
    n_iters : :obj:`int`, optional
        Number of iterations for each Markov Chain after burnin, by default 500
    n_burnin : :obj:`int`, optional
        Number of iterations for burnin of each Markov Chain, by default 100
    return_tree : :obj:`bool`, optional
        Return the inferred cell-lineage tree, by default False
    experiment : :obj:`bool`, optional
        Is in the experiment mode (the log won't be shown), by default False
    jar_path : :obj:`str`, optional
        Explicit path to ``SiCloneFiTComplete.jar``. The JAR only needs to be
        readable; it does not need executable permissions.
    java_executable : :obj:`str`, optional
        Explicit path to the Java runtime, by default ``java`` from the configured
        external tools directories or ``PATH``.
    tools_dir : :obj:`str`, optional
        Directory to search before ``SCPHYLO_TOOLS_DIR`` and ``PATH``.

    Returns
    -------
    :class:`pandas.DataFrame`
        A conflict-free matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1) and absence (0).
    """
    jar_path = scp.ul.resolve_external_file(
        jar_path or "SiCloneFiTComplete.jar",
        "SiCloneFit",
        tools_dir=tools_dir,
        path_option="jar_path",
    )
    java_executable = scp.ul.resolve_executable(
        java_executable or "java",
        "SiCloneFit",
        tools_dir=tools_dir,
        path_option="java_executable",
    )

    if not experiment:
        scp.logg.info(
            f"running SiCloneFit with alpha={alpha}, beta={beta}, n_iters={n_iters}"
        )

    with scp.ul.tmpdirsys(suffix=".siclonefit") as tmpdirname:
        workdir = Path(tmpdirname)
        input_path = workdir / "siclonefit.input"
        cellnames_path = workdir / "siclonefit.cellnames"
        genenames_path = workdir / "siclonefit.genenames"
        log_path = workdir / "siclonefit.log"

        df_input.T.reset_index(drop=True).to_csv(input_path, sep=" ", header=None)
        cellnames_path.write_text(" ".join(map(str, df_input.index)))
        genenames_path.write_text(" ".join(map(str, df_input.columns)))
        I_mtr = df_input.values

        command = [
            java_executable,
            "-jar",
            jar_path,
            "-m",
            str(df_input.shape[0]),
            "-n",
            str(df_input.shape[1]),
            "-ipMat",
            input_path,
            "-fp",
            str(alpha),
            "-fn",
            str(beta),
            "-df",
            "0",
            "-missing",
            str(np.sum(I_mtr == 3) / I_mtr.size),
            "-iter",
            str(n_iters),
            "-cellNames",
            cellnames_path,
            "-geneNames",
            genenames_path,
            "-r",
            str(n_restarts),
            "-burnin",
            str(n_burnin),
            "-outDir",
            workdir,
        ]
        s_time = time.perf_counter()
        with log_path.open("w") as log:
            scp.ul.run_external(
                command,
                "SiCloneFit",
                stdout=log,
                stderr=subprocess.STDOUT,
                log_path=log_path,
            )
        running_time = time.perf_counter() - s_time

        out_dirs = sorted(workdir.glob("*samples/best"))
        if not out_dirs:
            raise scp.ul.ExternalToolExecutionError(
                "SiCloneFit completed without producing its `samples/best` output."
            )
        out_dir = out_dirs[0]
        genotype_path = out_dir / "best_MAP_predicted_genotype.txt"
        try:
            df_output = pd.read_csv(
                genotype_path,
                sep=r"\s+",
                header=None,
                index_col=0,
            ).T
        except (OSError, pd.errors.ParserError) as error:
            raise scp.ul.ExternalToolExecutionError(
                "SiCloneFit produced an unreadable genotype matrix."
            ) from error
        if df_output.shape != df_input.shape:
            raise scp.ul.ExternalToolExecutionError(
                "SiCloneFit produced a genotype matrix with shape "
                f"{df_output.shape}; expected {df_input.shape}."
            )
        df_output.columns = df_input.columns.copy()
        df_output.index = df_input.index.copy()
        df_output.index.name = "cellIDxmutID"

        tree = None
        if return_tree and not experiment:
            tree_path = out_dir / "best_MAP_tree.txt"
            try:
                tree = tree_path.read_text().splitlines()[0].strip()
            except (OSError, IndexError) as error:
                raise scp.ul.ExternalToolExecutionError(
                    "SiCloneFit did not produce a readable MAP tree."
                ) from error

    if not experiment:
        scp.ul.stat(df_input, df_output, alpha, beta, running_time)
        if return_tree:
            return df_output, tree
        else:
            return df_output
    else:
        is_cf = scp.ul.is_conflict_free_gusfield(df_output)
        nll = scp.ul.calc_nll_matrix(df_input, df_output, alpha, beta)
        return df_output, running_time, is_cf, nll
