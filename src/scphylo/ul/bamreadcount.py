#!/usr/bin/env python
import glob
import sys

import anndata as ad
import pandas as pd


def get_nucleotides(line):
    ref = line.split()[2].upper()
    char = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
    char[ref] = int(line.split()[4].split(":")[1])
    char[line.split()[5].split(":")[0]] = int(line.split()[5].split(":")[1])
    char[line.split()[6].split(":")[0]] = int(line.split()[6].split(":")[1])
    char[line.split()[7].split(":")[0]] = int(line.split()[7].split(":")[1])
    char[line.split()[8].split(":")[0]] = int(line.split()[8].split(":")[1])
    cov = sum(char[a] for a in ["A", "C", "G", "T", "N"])
    return ref, char, cov


def is_mut(x):
    g = 3
    if x["COV"] > 0:
        g = 0
    if x[x["Allele"][0]] > 0:
        g = 1
    return g, x[x["Allele"][0]]


outdir = sys.argv[1]
vcf = pd.read_csv(sys.argv[2])
data = []
for file in glob.glob(f"{outdir}/*.bamreadcount"):
    with open(file) as fin:
        name = file.split("/")[-1].replace(".bamreadcount", "")
        for line in fin:
            line = line.strip()
            a = ".".join(line.split("\t")[:2])
            ref, char, cov = get_nucleotides(line)
            res = {"CELL": name, "MUT": a, "REF": ref, "COV": cov}
            for k, v in char.items():
                res[k] = v
            data.append(res)
df = pd.DataFrame(data)

vcf["MUT"] = vcf["CHROM"] + "." + vcf["POS"].astype(str)
df = pd.merge(df, vcf, how="left", left_on="MUT", right_on="MUT")
df["IS_MUT"] = df.apply(lambda x: is_mut(x)[0], axis=1).astype(int)
df["MUTANT"] = df.apply(lambda x: is_mut(x)[1], axis=1).astype(int)

sc = df.pivot_table(index="CELL", columns="MUT", values="IS_MUT")
adata = ad.AnnData(sc)
adata.layers["mutant"] = df.pivot_table(index="CELL", columns="MUT", values="MUTANT")
adata.layers["total"] = df.pivot_table(index="CELL", columns="MUT", values="COV")
adata.var = pd.merge(adata.var, vcf, how="left", left_index=True, right_on="MUT")
adata.write(outdir + "/_bamreadcount.h5ad.gz", compression="gzip")
