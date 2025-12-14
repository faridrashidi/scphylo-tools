#!/usr/bin/env python
import glob
import os
import sys

import anndata as ad
import pandas as pd

outdir = sys.argv[1]
adata = ad.AnnData()
files = glob.glob(outdir + "/*.transcript.bam")
kinds = {"expected_count": "count", "TPM": "tpm", "FPKM": "fpkm"}
for kind in kinds.keys():
    i = 0
    expr = pd.DataFrame()
    for file in files:
        file = file[: -len(".transcript.bam")] + ".genes.results"
        df = pd.read_csv(file, sep="\t")
        name = os.path.basename(file)[: -len(".genes.results")]
        df = df[["gene_id", "expected_count", "TPM", "FPKM"]]
        df = df.set_index(["gene_id"])
        if i == 0:
            expr = pd.DataFrame(index=df.index)
        expr[name] = df[kind]
        i += 1
    expr = expr.transpose()
    expr = expr.sort_index(ascending=True)
    if kind == "expected_count":
        adata = ad.AnnData(expr)
    else:
        adata.layers[kind] = expr.values
adata.write(outdir + "/_expression.h5ad.gz", compression="gzip")
