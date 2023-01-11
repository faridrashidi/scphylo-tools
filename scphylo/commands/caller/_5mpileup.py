import glob
import os

import click
import pandas as pd

import scphylo as scp
from scphylo.ul._servers import cmd, write_cmds_get_main


@click.command(short_help="Run mpileup.")
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
    default="5-00:00:00",
    type=str,
    show_default=True,
    help="Time.",
)
@click.option(
    "--mem",
    default="150",
    type=str,
    show_default=True,
    help="Memory.",
)
@click.option(
    "--afterok",
    default=None,
    type=str,
    show_default=True,
    help="Afterok.",
)
def mpileup(outdir, ref, time, mem, afterok):
    """Run mpileup.

    scphylo caller mpileup /path/to/in/dir hg19|hg38|mm10

    BAM files (*.markdup_bqsr.bam) --> VCF file (_pileupcalls.vcf)
    """
    if ref == "hg19":
        config = scp.settings.hg19
    elif ref == "mm10":
        config = scp.settings.mm10
    elif ref == "hg38":
        config = scp.settings.hg38
    elif ref == "m10x":
        config = scp.settings.m10x
    else:
        config = None

    def cmds(samples):
        cmds = ""
        cmds += cmd([f"module load {scp.settings.tools['bcftools']}"])
        files = " ".join([f"{outdir}/{s}.markdup_bqsr.bam" for s in samples])
        cmds += cmd(
            [
                "bcftools mpileup",
                "--min-BQ 30",
                "--min-MQ 0",
                "--count-orphans",
                "--ignore-overlaps",
                "--output-type u",
                "--max-depth 5000",
                '--annotate "DP,AD"',
                "--threads $SLURM_CPUS_PER_TASK",
                f"-f {config['ref']}",
                files,
                "|",
                "bcftools call",
                "--output-type v",
                # "---multiallelic-caller",
                "--consensus-caller",
                "--variants-only",
                f"-o {outdir}/_pileupcalls.vcf",
            ]
        )
        cmds += cmd(["echo Done!"], islast=True)
        return cmds

    files = glob.glob(f"{outdir}/*.markdup_bqsr.bam")
    temp = []
    for file in files:
        file = file.split("/")[-1].replace(".markdup_bqsr.bam", "")
        temp.append(file)
    df_cmds = pd.DataFrame()
    df_cmds["cmd"] = [cmds(temp)]

    cmdmain = write_cmds_get_main(
        df_cmds,
        "mpileup",
        time,
        mem,
        None,
        16,
        scp.settings.tools["email"],
        f"{outdir}/_tmp",
        afterok,
    )
    os.system(cmdmain)

    return None
