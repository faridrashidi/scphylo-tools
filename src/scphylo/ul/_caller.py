import glob

import pandas as pd

import scphylo as scp


def is_paired(indir):
    files = glob.glob(f"{indir}/*.fastq.gz")
    last_char_set = {x[-len(".fastq.gz") - 2 : -len(".fastq.gz")] for x in files}
    temp = []
    if last_char_set == {"_1", "_2"}:
        for file in files:
            file = file.split("/")[-1].replace(".fastq.gz", "")
            if file[-2:] == "_1":
                temp.append({"sample": file[:-2], "is_paired": True})
    else:
        for file in files:
            file = file.split("/")[-1].replace(".fastq.gz", "")
            temp.append({"sample": file, "is_paired": False})
    return temp


def get_samples_df(outdir, normal=None):
    files = glob.glob(f"{outdir}/*.markdup_bqsr.bam")
    temp = []
    for file in files:
        file = file.split("/")[-1].replace(".markdup_bqsr.bam", "")
        if normal is not None:
            if file != normal:
                temp.append({"sample": file})
        else:
            temp.append({"sample": file})
    df_cmds_tmp = pd.DataFrame(temp)
    if normal is not None:
        scp.logg.info(
            f"Tumor samples: {','.join(df_cmds_tmp['sample'].values)} & Normal sample:"
            f" {normal}"
        )
    return df_cmds_tmp
