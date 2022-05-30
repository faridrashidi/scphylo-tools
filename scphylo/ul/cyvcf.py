#!/usr/bin/env python
import sys

import anndata as ad
import numpy as np
import pandas as pd
from cyvcf2 import VCF

import scphylo as scp

infile = sys.argv[1]
vcf = VCF(infile)
info = vcf.get_header_type("ANN")["Description"].split(" | ")
info[0] = "Allele"
info[-1] = "ERRORS / WARNINGS / INFO"
info = ["CHROM", "POS", "REF", "ALT", "START", "END"] + info

cells = vcf.samples
muts = []
m_gen = []
m_ref = []
m_alt = []
m_cna = []
for var in VCF(infile):
    if var.FILTER is None:
        if var.is_snp or var.is_indel:
            row = [var.CHROM, var.POS, var.REF, var.ALT, var.start + 1, var.end]
            if var.INFO.get("ANN") is not None:
                ann = var.INFO.get("ANN").split(",")[0].split("|")
                row += ann
            m_gen.append(var.gt_types.tolist())
            m_ref.append(var.gt_ref_depths.tolist())
            m_alt.append(var.gt_alt_depths.tolist())
            muts.append(row)
            cnas_per_var = []

m_gen = np.array(m_gen)
m_ref = np.array(m_ref)
m_alt = np.array(m_alt)
muts = pd.DataFrame(muts, columns=info)
muts.index = muts.index.map(lambda x: f"mut{x}")
cells = pd.DataFrame(index=cells)

adata = ad.AnnData(np.zeros((len(muts), len(cells))))
adata.obs = muts
adata.var = cells
if m_gen.shape[1] > 0:
    adata.layers["genotype"] = m_gen
adata.layers["total"] = m_ref + m_alt
adata.layers["mutant"] = m_alt
adata = adata.T
outdir, basename = scp.ul.dir_base(infile)

bad = [
    "Annotation_Impact",
    "Gene_ID",
    "Feature_ID",
    "Rank",
    "cDNA.pos / cDNA.length",
    "CDS.pos / CDS.length",
    "AA.pos / AA.length",
    "Distance",
    "ERRORS / WARNINGS / INFO",
]
adata.var.drop(bad, axis=1, inplace=True)
if len(adata.obs_names) == 2:
    for obs in adata.obs_names:
        adata.var[f"{obs}_TOTAL_READS"] = adata[obs].layers["total"][0]
        adata.var[f"{obs}_MUTANT_READS"] = adata[obs].layers["mutant"][0]
adata.var.to_csv(outdir + f"/{basename}.tsv", sep="\t")
adata.write(outdir + f"/_{basename[:-len('.ann')]}.h5ad.gz", compression="gzip")
