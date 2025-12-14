import glob
import os

import click
import pandas as pd

import scphylo as scp
from scphylo.ul._servers import cmd, write_cmds_get_main


@click.command(short_help="Run deFuse.")
@click.argument(
    "indir",
    required=True,
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, readable=True, resolve_path=True
    ),
)
@click.argument(
    "outdir",
    required=True,
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, readable=True, resolve_path=True
    ),
)
@click.argument(
    "ref",
    required=True,
    type=click.Choice(scp.settings.refs),
)
@click.option(
    "--time",
    default="0-02:00:00",
    type=str,
    show_default=True,
    help="Time.",
)
@click.option(
    "--mem",
    default="50",
    type=str,
    show_default=True,
    help="Memory.",
)
@click.option(
    "--n_threads",
    default="16",
    type=str,
    show_default=True,
    help="Number of Threads.",
)
@click.option(
    "--afterok",
    default=None,
    type=str,
    show_default=True,
    help="Afterok.",
)
def defuse(indir, outdir, ref, time, mem, n_threads, afterok):
    """Run deFuse.

    scphylo caller defuse ...

    FastQ files (*.fastq) -> deFuse files
    """
    if ref == "hg19":
        config = scp.settings.hg19
    elif ref == "mm10":
        config = scp.settings.mm10
    elif ref == "hg38":
        config = scp.settings.hg38
    else:
        config = None

    def cmds(sample):
        cmds = ""
        cmds += cmd([f"module load {scp.settings.tools['defuse']}"])
        cmds += cmd([f"mkdir -p {outdir}/{sample}"])
        cmds += cmd(
            [
                "defuse.pl",
                f"-c {config['defuse_config']}",
                f"-o {outdir}/{sample}",
                f"-d {config['defuse_db']}",
                f"-1 {indir}/{sample}_R1_001.fastq",
                f"-2 {indir}/{sample}_R2_001.fastq",
                "-s direct",
                "-p $SLURM_CPUS_PER_TASK",
            ]
        )
        cmds += cmd(["echo Done!"], islast=True)
        return cmds

    df = pd.DataFrame()
    temp = []
    for x in glob.glob(f"{indir}/*_R1_001.fastq"):
        temp.append(x.split("/")[-1].replace("_R1_001.fastq", ""))
    df["sample"] = temp
    df["cmd"] = df.apply(lambda x: cmds(x["sample"]), axis=1)
    scp.logg.info(f"Number of samples: {df.shape[0]}")

    cmdmain = write_cmds_get_main(
        df,
        "defuse",
        time,
        mem,
        None,
        n_threads,
        scp.settings.tools["email"],
        f"{outdir}/_tmp",
        afterok,
    )
    os.system(cmdmain)

    return None
