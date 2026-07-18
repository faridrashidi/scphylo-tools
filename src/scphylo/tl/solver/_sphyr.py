import time
from pathlib import Path

import pandas as pd

import scphylo as scp


def sphyr(
    df_input,
    alpha,
    beta,
    n_restarts=10,
    n_threads=1,
    time_limit=None,
    n_cell_clusters=10,
    n_mut_clusters=15,
    executable=None,
    tools_dir=None,
):
    """Solving using SPhyR.

    Tumor phylogeny estimation from single-cell sequencing data under loss and error
    :cite:`SPhyR`.

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
        Number of threads, by default 10
    n_threads : :obj:`int`, optional
        Number of threads, by default 1
    time_limit : :obj:`int`, optional
        Time limit (in seconds), by default None
    n_cell_clusters : :obj:`int`, optional
        Number of cell clusters, by default 10
    n_mut_clusters : :obj:`int`, optional
        Number of mutation clusters, by default 15
    executable : :obj:`str`, optional
        Explicit path to ``sphyr_kDPFC``. By default, search the configured external
        tools directories and ``PATH``.
    tools_dir : :obj:`str`, optional
        Directory to search before ``SCPHYLO_TOOLS_DIR`` and ``PATH``.

    Returns
    -------
    :class:`pandas.DataFrame`
        A conflict-free matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1) and absence (0).
    """
    executable = scp.ul.resolve_executable(
        executable or "sphyr_kDPFC",
        "SPhyR",
        tools_dir=tools_dir,
        path_option="executable",
    )

    scp.logg.info(
        f"running SPhyR with alpha={alpha}, beta={beta}, n_restarts={n_restarts}, "
        f"n_threads={n_threads}, time_limit={time_limit}, "
        f"n_cell_clusters={n_cell_clusters}, n_mut_clusters={n_mut_clusters}"
    )

    with scp.ul.tmpdirsys(suffix=".sphyr") as tmpdirname:
        workdir = Path(tmpdirname)
        input_path = workdir / "sphyr.input"
        output_path = workdir / "sphyr.output"
        log_path = workdir / "sphyr.log"
        cellnames_path = workdir / "sphyr.cellnames"
        mutnames_path = workdir / "sphyr.mutnames"

        with input_path.open("w") as fout:
            fout.write(f"{df_input.shape[0]} #cells\n{df_input.shape[1]} #SNVs\n")
            df_input.replace(3, -1).to_csv(fout, sep=" ", header=None, index=None)
        cellnames_path.write_text("\n".join(map(str, df_input.index)) + "\n")
        mutnames_path.write_text("\n".join(map(str, df_input.columns)) + "\n")

        command = [
            executable,
            input_path,
            "-a",
            str(alpha),
            "-b",
            str(beta),
            "-N",
            str(n_restarts),
            "-t",
            str(n_threads),
            "-lC",
            str(n_mut_clusters),
            "-lT",
            str(n_cell_clusters),
            "-T",
            str(time_limit if time_limit is not None else -1),
            "-k",
            "0",
        ]

        s_time = time.perf_counter()
        with output_path.open("w") as output, log_path.open("w") as log:
            scp.ul.run_external(
                command,
                "SPhyR",
                stdout=output,
                stderr=log,
                log_path=log_path,
            )
        running_time = time.perf_counter() - s_time

        if output_path.stat().st_size == 0:
            raise scp.ul.ExternalToolExecutionError(
                "SPhyR completed without producing a genotype matrix."
            )
        try:
            df_output = pd.read_csv(
                output_path,
                sep=r"\s+",
                skiprows=[0, 1],
                header=None,
            )
        except (OSError, pd.errors.ParserError) as error:
            raise scp.ul.ExternalToolExecutionError(
                "SPhyR produced an unreadable genotype matrix."
            ) from error
        if df_output.shape != df_input.shape:
            raise scp.ul.ExternalToolExecutionError(
                "SPhyR produced a genotype matrix with shape "
                f"{df_output.shape}; expected {df_input.shape}."
            )
        df_output.index = df_input.index.copy()
        df_output.columns = df_input.columns.copy()
        df_output.index.name = "cellIDxmutID"

    scp.ul.stat(df_input, df_output, alpha, beta, running_time)

    return df_output
