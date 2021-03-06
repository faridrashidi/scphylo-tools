#!/usr/bin/env python
import glob
import sys

import pandas as pd

import scphylo as scp


def get_numreads_percmapped(file):
    data = {
        "total_reads": None,
        "uniquely_mapped": None,
        "uniquely_mapped_percent": None,
        "num_splices": None,
        "num_GCAG_splices": None,
        "insertion_length": None,
        "deletion_length": None,
        "unmapped_tooshort_percent": None,
        "avg_mapped_read_length": None,
        "deletion_rate": None,
        "mismatch_rate": None,
        "avg_input_read_length": None,
        "num_ATAC_splices": None,
        "num_annotated_splices": None,
        "num_GTAG_splices": None,
        "multimapped_toomany": None,
        "unmapped_mismatches": None,
        "unmapped_mismatches_percent": None,
        "unmapped_other": None,
        "insertion_rate": None,
        "unmapped_other_percent": None,
        "multimapped_percent": None,
        "multimapped": None,
        "num_noncanonical_splices": None,
        "unmapped_tooshort": None,
        "multimapped_toomany_percent": None,
    }
    with open(file) as fin:
        for line in fin:
            if "Uniquely mapped reads % |" in line:
                data["uniquely_mapped_percent"] = float(
                    line.replace("Uniquely mapped reads % |", "")
                    .strip()
                    .replace("%", "")
                )
            if "Number of splices: Total |" in line:
                data["num_splices"] = int(
                    line.replace("Number of splices: Total |", "").strip()
                )
            if "Number of splices: GC/AG |" in line:
                data["num_GCAG_splices"] = int(
                    line.replace("Number of splices: GC/AG |", "").strip()
                )
            if "Insertion average length |" in line:
                data["insertion_length"] = float(
                    line.replace("Insertion average length |", "")
                    .strip()
                    .replace("%", "")
                )
            if "Deletion average length |" in line:
                data["deletion_length"] = float(
                    line.replace("Deletion average length |", "")
                    .strip()
                    .replace("%", "")
                )
            if "% of reads unmapped: too short |" in line:
                data["unmapped_tooshort_percent"] = float(
                    line.replace("% of reads unmapped: too short |", "")
                    .strip()
                    .replace("%", "")
                )
            if "Average mapped length |" in line:
                data["avg_mapped_read_length"] = float(
                    line.replace("Average mapped length |", "").strip().replace("%", "")
                )
            if "Deletion rate per base |" in line:
                data["deletion_rate"] = float(
                    line.replace("Deletion rate per base |", "")
                    .strip()
                    .replace("%", "")
                )
            if "Mismatch rate per base, % |" in line:
                data["mismatch_rate"] = float(
                    line.replace("Mismatch rate per base, % |", "")
                    .strip()
                    .replace("%", "")
                )
            if "Average input read length |" in line:
                data["avg_input_read_length"] = float(
                    line.replace("Average input read length |", "")
                    .strip()
                    .replace("%", "")
                )
            if "Number of splices: AT/AC |" in line:
                data["num_ATAC_splices"] = int(
                    line.replace("Number of splices: AT/AC |", "").strip()
                )
            if "Number of splices: Annotated (sjdb) |" in line:
                data["num_annotated_splices"] = int(
                    line.replace("Number of splices: Annotated (sjdb) |", "").strip()
                )
            if "Number of splices: GT/AG |" in line:
                data["num_GTAG_splices"] = int(
                    line.replace("Number of splices: GT/AG |", "").strip()
                )
            if "Uniquely mapped reads number |" in line:
                data["uniquely_mapped"] = int(
                    line.replace("Uniquely mapped reads number |", "").strip()
                )
            if "Number of reads mapped to too many loci |" in line:
                data["multimapped_toomany"] = int(
                    line.replace(
                        "Number of reads mapped to too many loci |", ""
                    ).strip()
                )
            if "Number of reads unmapped: too many mismatches |" in line:
                data["unmapped_mismatches"] = int(
                    line.replace(
                        "Number of reads unmapped: too many mismatches |", ""
                    ).strip()
                )
            if "% of reads unmapped: too many mismatches |" in line:
                data["unmapped_mismatches_percent"] = float(
                    line.replace("% of reads unmapped: too many mismatches |", "")
                    .strip()
                    .replace("%", "")
                )
            if "Number of input reads |" in line:
                data["total_reads"] = int(
                    line.replace("Number of input reads |", "").strip()
                )
            if "Number of reads unmapped: other |" in line:
                data["unmapped_other"] = int(
                    line.replace("Number of reads unmapped: other |", "").strip()
                )
            if "Insertion rate per base |" in line:
                data["insertion_rate"] = float(
                    line.replace("Insertion rate per base |", "")
                    .strip()
                    .replace("%", "")
                )
            if "% of reads unmapped: other |" in line:
                data["unmapped_other_percent"] = float(
                    line.replace("% of reads unmapped: other |", "")
                    .strip()
                    .replace("%", "")
                )
            if "% of reads mapped to multiple loci |" in line:
                data["multimapped_percent"] = float(
                    line.replace("% of reads mapped to multiple loci |", "")
                    .strip()
                    .replace("%", "")
                )
            if "Number of reads mapped to multiple loci |" in line:
                data["multimapped"] = int(
                    line.replace(
                        "Number of reads mapped to multiple loci |", ""
                    ).strip()
                )
            if "Number of splices: Non-canonical |" in line:
                data["num_noncanonical_splices"] = int(
                    line.replace("Number of splices: Non-canonical |", "").strip()
                )
            if "Number of reads unmapped: too short |" in line:
                data["unmapped_tooshort"] = float(
                    line.replace("Number of reads unmapped: too short |", "")
                    .strip()
                    .replace("%", "")
                )
            if "% of reads mapped to too many loci |" in line:
                data["multimapped_toomany_percent"] = float(
                    line.replace("% of reads mapped to too many loci |", "")
                    .strip()
                    .replace("%", "")
                )
    return data


outdir = sys.argv[1]
files = glob.glob(outdir + "/*.star.log")
data = []
samples = []
for file in files:
    _, sample = scp.ul.dir_base(file)
    samples.append(sample[: -len(".star")])
    d = get_numreads_percmapped(f"{outdir}/{sample}.log")
    data.append(d)
df = pd.DataFrame(data, index=samples)
df.to_csv(outdir + "/_starinfo.csv")
