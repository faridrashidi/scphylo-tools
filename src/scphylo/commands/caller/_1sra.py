import glob
import os
import subprocess

import click
import pandas as pd

import scphylo as scp
from scphylo.ul._servers import cmd, write_cmds_get_main


@click.command(short_help="Run SRA.")
@click.argument(
    "srafile",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument(
    "outdir",
    required=True,
    type=click.Path(file_okay=False, dir_okay=True, readable=True, resolve_path=True),
)
@click.option(
    "--check",
    is_flag=True,
    default=False,
    type=bool,
    show_default=True,
    help="Check the unfinished jobs on biowulf.",
)
def sra(srafile, outdir, check):
    """Run SRA.

    scphylo caller sra path/to/SraRunTable.csv path/to/output/directory

    Columns name of the csv file: `Library Name`, `Run`, `LibraryLayout`
    """
    scp.ul.mkdir(outdir)

    def cmds(srr_id, name, layout):
        cmds = ""
        if layout == "PAIRED":
            cmds += cmd(
                [
                    "fastq-dump",
                    "--split-files",
                    f"--outdir {outdir}",
                    f"{srr_id}",
                ]
            )
            if srr_id != name:
                cmds += cmd(
                    [
                        "mv",
                        f"{outdir}/{srr_id}_1.fastq",
                        f"{outdir}/{name}_1.fastq",
                    ]
                )
                cmds += cmd(
                    [
                        "mv",
                        f"{outdir}/{srr_id}_2.fastq",
                        f"{outdir}/{name}_2.fastq",
                    ]
                )
            cmds += cmd(
                [
                    "gzip",
                    f"{outdir}/{name}_1.fastq",
                ]
            )
            cmds += cmd(
                [
                    "gzip",
                    f"{outdir}/{name}_2.fastq",
                ]
            )
        elif layout == "SINGLE":
            cmds += cmd(
                [
                    "fastq-dump",
                    f"--outdir {outdir}",
                    f"{srr_id}",
                ]
            )
            if srr_id != name:
                cmds += cmd(
                    [
                        "mv",
                        f"{outdir}/{srr_id}.fastq",
                        f"{outdir}/{name}.fastq",
                    ]
                )
            cmds += cmd(
                [
                    "gzip",
                    f"{outdir}/{name}.fastq",
                ]
            )

        cmds += cmd(["echo Done!"], islast=True)
        return cmds

    if not check:
        df = pd.read_csv(srafile, low_memory=False)
        scp.logg.info(f"There are {df.shape[0]} samples.")

        if df["Library Name"].nunique() != df.shape[0]:
            scp.logg.error(
                "Number of unique `Library Name` and total samples are not equal!"
            )

        df["cmd"] = df.apply(
            lambda x: cmds(x["Run"], x["Library Name"], x["LibraryLayout"]), axis=1
        )
        cmdmain = write_cmds_get_main(
            df,
            "sra",
            "12:00:00",
            "10",
            "python/3.7,sratoolkit/2.10.8",
            1,
            "farid.rsh@gmail.com",
            f"{outdir}/_tmp",
            None,
        )
        os.system(cmdmain)
    else:
        n = 0
        for file in glob.glob(f"{outdir}/_tmp/log/_sra/*.o"):
            out = subprocess.getoutput(f"cat {file}")
            if out.count("\nDone!") != out.count("(") and out.count(
                "\nDone!"
            ) != out.count(")"):
                with open(file) as fin:
                    scp.logg.info(
                        fin.readlines()[1].replace(")", "").replace("(", "").strip()
                    )
                    scp.logg.info("\n")
                    n += 1
        scp.logg.info(f"There are {n} unfinished jobs.")

    return None
