import os

import anndata as ad
import ete3
import networkx as nx
import numpy as np
import pandas as pd

import scphylo as scp


def read(filepath):
    """Read genotype matrix and read-count matrix.

    The genotype matrix must be in the in format of :class:`pandas.DataFrame`
    The read-count matrix must be in the format of :class:`anndata.AnnData`.

    Parameters
    ----------
    filepath : :obj:`str`
        The path to the file. The extension must be one of
        [`.tsv`, `.SC`, `.CFMatrix`, `.h5ad`, `.h5ad.gz`, `.nwk`]

    Returns
    -------
    :class:`pandas.DataFrame` or :class:`anndata.AnnData`
        Depends on the format of the input file the output type is different.
    """
    ext = os.path.splitext(filepath)[-1]
    if ext in [".SC", ".CFMatrix", ".before_FP_FN_NA", ".tsv"]:
        sc = pd.read_table(filepath, index_col=0)
        if len(sc.columns) != len(set(sc.columns)):
            scp.logg.error("Mutation ids must be unique!")
        return sc
    elif ext in [".csv"]:
        sc = pd.read_csv(filepath, index_col=0)
        if len(sc.columns) != len(set(sc.columns)):
            scp.logg.error("Mutation ids must be unique!")
        return sc
    elif ext in [".h5ad", ".gz"]:
        return ad.read(filepath)
    elif ext in [".nwk"]:
        return _read_nwk(filepath)
    else:
        scp.logg.error("Extension is wrong!")
        return None


def write(obj, filepath):
    """Write genotype matrix or read-count matrix into a file.

    Parameters
    ----------
    obj : :class:`pandas.DataFrame` or :class:`anndata.AnnData`
        The input object which is going to be written in a file.
    filepath : :obj:`str`
        The file path where the `obj` must be written in.
    """
    if isinstance(obj, pd.DataFrame):
        obj.index.name = "cellIDxmutID"
        obj.to_csv(filepath, sep="\t")
    elif isinstance(obj, ad.AnnData):
        obj.write(filepath + ".h5ad.gz", compression="gzip")
    else:
        scp.logg.error("Object instance is wrong!")


def _read_nwk(filepath):
    tree = ete3.Tree(filepath, format=1)
    G = nx.DiGraph()
    node2id = {}
    i = 0
    for n in tree.traverse("postorder"):
        if n.name == "" or "Inner" in n.name:
            G.add_node(i, label="––")
        else:
            G.add_node(i, label=str(n.name))
        node2id[n] = i
        i += 1

    for p in tree.traverse("postorder"):
        pn = node2id[p]
        for c in p.children:
            cn = node2id[c]
            G.add_edge(pn, cn)

    i = 0
    for e, u, _ in G.edges.data("label"):
        G.edges[(e, u)]["label"] = f"m{i}"
        i += 1
    G.graph["normal_cells"] = []
    G.graph["splitter_mut"] = "\n"
    G.graph["splitter_cell"] = "\n"
    data = scp.ul.to_cfmatrix(G)
    return data


def to_vcf(adata, filepath):
    def _calculate_pl(ref_count, alt_count):
        if ref_count == 0 and alt_count == 0:
            return [0, 0, 0]
        elif ref_count == 0:
            return [0, 0, 255]
        elif alt_count == 0:
            return [255, 0, 0]
        else:
            total_count = ref_count + alt_count
            p_ref = ref_count / total_count
            p_alt = alt_count / total_count
            p_ref_given_data = p_ref
            p_alt_given_data = p_alt
            p_ref_ref_given_data = p_ref * p_ref
            pl_ref = -10 * np.math.log10(p_ref_given_data)
            pl_alt = -10 * np.math.log10(p_alt_given_data)
            pl_ref_ref = -10 * np.math.log10(p_ref_ref_given_data)
            return [int(pl_ref), int(pl_ref_ref), int(pl_alt)]

    gt_mtx = np.where(
        adata.X == 0,
        "0/0",
        np.where(adata.X == 1, "0/1", np.where(adata.X == 3, "./.", adata.X)),
    )
    ref_mtx = adata.layers["total"] - adata.layers["mutant"]
    alt_mtx = adata.layers["mutant"]
    with open(filepath, "w") as fout:
        fout.write("##fileformat=VCFv4.3\n")
        fout.write("##contig=<ID=chr1,length=1000000>\n")
        fout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fout.write(
            '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized,'
            " Phred-scaled likelihoods for genotypes as defined in the VCF"
            ' specification">\n'
        )
        fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        fout.write("\t".join(adata.obs_names) + "\n")
        for j in range(adata.shape[1]):
            samples = []
            for i in range(adata.shape[0]):
                pl = _calculate_pl(ref_mtx[i, j], alt_mtx[i, j])
                samples.append(gt_mtx[i, j] + ":" + ",".join(map(str, pl)))
            variant = (
                "\t".join(
                    [
                        "chr1",
                        str(j + 1),
                        ".",
                        "A",
                        "C",
                        ".",
                        "PASS",
                        ".",
                        "GT:PL",
                        "\t".join(samples),
                    ]
                )
                + "\n"
            )
            fout.write(variant)
