import os
import time

import scphylo as scp


def cellphy(adata, mode="fast"):
    """Solving using CellPhy.

    Accurate and fast probabilistic inference of single-cell phylogenies
    from scDNA-seq data :cite:`CellPhy`.

    Parameters
    ----------
    adata : :class:`anndata.AnnData`
        Input data contains layers of mutant and total.
    mode : :obj:`str`
        CellPhy mode, possible values are:

            - `fast`
            - `full`

    Returns
    -------
    :class:`pandas.DataFrame`
        A conflict-free matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1) and absence (0).
    """
    executable = scp.ul.executable("cellphy.sh", "CellPhy")

    if adata.shape[1] < 37:
        scp.logg.error("number of mutations must be greater than 37!")

    scp.logg.info(f"running CellPhy on {mode} mode")

    tmpdir = scp.ul.tmpdirsys(suffix=".cellphy")

    scp.io.to_vcf(adata, filepath=f"{tmpdir.name}/cellphy.vcf")

    cmd = f"{executable} {mode.upper()} {tmpdir.name}/cellphy.vcf"

    s_time = time.time()
    os.system(cmd)
    e_time = time.time()
    running_time = e_time - s_time
    type(running_time)

    with open(f"{tmpdir.name}/cellphy.vcf.raxml.bestTree") as fin:
        tree = fin.readline().strip()

    tmpdir.cleanup()

    return tree
