import subprocess

import click
import pandas as pd

import scphylo as scp
from scphylo.ul._servers import cmd, write_cmds_get_main


@click.command(short_help="Run fastp.")
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
    type=click.Path(file_okay=False, dir_okay=True, readable=True, resolve_path=True),
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
    default="10",
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
def fastp(indir, outdir, time, mem, afterok):
    """Run fastp.

    scphylo caller fastp /path/to/in/dir /path/to/out/dir

    FastQ files --> FastQ files
    """
    scp.ul.mkdir(outdir)

    def get_command(sample, is_paired):
        cmds = ""
        cmds += cmd([f"module load {scp.settings.tools['fastp']}"])
        if is_paired:
            infqs = [
                f"-i {indir}/{sample}_1.fastq.gz",
                f"-I {indir}/{sample}_2.fastq.gz",
            ]
            outfqs = [
                f"-o {outdir}/{sample}_1.fastq.gz",
                f"-O {outdir}/{sample}_2.fastq.gz",
            ]
        else:
            infqs = [
                f"-i {indir}/{sample}.fastq.gz",
            ]
            outfqs = [
                f"-o {outdir}/{sample}.fastq.gz",
            ]
        cmds += cmd(
            [
                "fastp",
                "-w 8",
                f"{' '.join(infqs)}",
                f"{' '.join(outfqs)}",
                f"-h {outdir}/{sample}.fastp.html",
                f"-j {outdir}/{sample}.fastp.json",
                f"> {outdir}/{sample}.fastp.log",
            ]
        )
        cmds += cmd(
            [
                "rm -rf",
                f"{outdir}/{sample}.fastp.html",
                f"{outdir}/{sample}.fastp.json",
                f"{outdir}/{sample}.fastp.log",
            ]
        )
        cmds += cmd(["echo Done!"], islast=True)
        return cmds

    temp = scp.ul.is_paired(indir)
    df_cmds = pd.DataFrame(temp)
    if df_cmds.shape[0] != df_cmds["sample"].nunique():
        scp.logg.error("Samples are not unique!")
    df_cmds["cmd"] = df_cmds.apply(
        lambda x: get_command(x["sample"], x["is_paired"]), axis=1
    )

    cmdmain = write_cmds_get_main(
        df_cmds,
        "fastp",
        time,
        mem,
        None,
        8,
        scp.settings.tools["email"],
        f"{outdir}/_tmp",
        afterok,
    )
    code = subprocess.getoutput(cmdmain)
    scp.logg.info(code)

    return None
