import subprocess

import click
import pandas as pd

import scphylo as scp
from scphylo.ul._servers import cmd, write_cmds_get_main


@click.command(short_help="Run BWA.")
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
@click.argument(
    "ref",
    required=True,
    type=click.Choice(scp.settings.refs),
)
@click.option(
    "--time",
    default="0-10:00:00",
    type=str,
    show_default=True,
    help="Time.",
)
@click.option(
    "--mem",
    default="90",
    type=str,
    show_default=True,
    help="Memory (mem must be at least `n_threads * 10GB`).",
)
@click.option(
    "--n_threads",
    default="8",
    type=str,
    show_default=True,
    help="Number of threads.",
)
@click.option(
    "--afterok",
    default=None,
    type=str,
    show_default=True,
    help="Afterok.",
)
@click.option(
    "--is_pdx",
    default=False,
    type=bool,
    show_default=True,
    is_flag=True,
    help="Is the the PDX.",
)
def bwa(indir, outdir, ref, time, mem, n_threads, afterok, is_pdx):
    """Run bwa-mem.

    scphylo caller bwa /path/to/in/dir /path/to/out/dir hg19|hg38|mm10 --is_pdx

    FastQ files (*.fastq.gz) --> BAM files (*.mapped.bam)

    if --is_pdx: FastQ files (*.fastq.gz) --> FastQ files (*.fastq.gz)
    """
    scp.ul.mkdir(outdir)
    if ref == "hg19":
        config = scp.settings.hg19
    elif ref == "mm10":
        config = scp.settings.mm10
    elif ref == "hg38":
        config = scp.settings.hg38
    else:
        config = None

    def get_command(sample, is_paired):
        cmds = ""
        cmds += cmd([f"module load {scp.settings.tools['bwa']}"])
        cmds += cmd([f"module load {scp.settings.tools['samtools']}"])
        if is_paired:
            infqs = [f"{indir}/{sample}_{i}.fastq.gz" for i in range(1, 3)]
        else:
            infqs = [f"{indir}/{sample}.fastq.gz"]
        cmds += cmd(
            [
                "bwa mem",
                f"-t {n_threads}",
                "-M",
                f"{config['bwa_ref']}",
                f"{' '.join(infqs)}",
                "|",
                f"samtools sort -@{n_threads} -m 10000000000",
                f"-o {outdir}/{sample}.mapped.bam -",
            ]
        )
        if is_pdx:
            cmds += cmd(
                [
                    "samtools view",
                    "-b -f 4",
                    f"{outdir}/{sample}.mapped.bam",
                    f"> {outdir}/{sample}.unmapped.bam",
                ]
            )
            cmds += cmd(
                [
                    "samtools fastq",
                    f"-1 {outdir}/{sample}_1.fastq",
                    f"-2 {outdir}/{sample}_2.fastq",
                    "-0 /dev/null -s /dev/null -n",
                    f"{outdir}/{sample}.unmapped.bam",
                ]
            )
            cmds += cmd(
                [
                    "gzip",
                    f"{outdir}/{sample}_1.fastq",
                ]
            )
            cmds += cmd(
                [
                    "gzip",
                    f"{outdir}/{sample}_2.fastq",
                ]
            )
            cmds += cmd(
                [
                    "rm -rf",
                    f"{outdir}/{sample}.mapped.bam",
                    f"{outdir}/{sample}.unmapped.bam",
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
        "bwa",
        time,
        mem,
        None,
        n_threads,
        scp.settings.tools["email"],
        f"{outdir}/_tmp",
        afterok,
    )
    code = subprocess.getoutput(cmdmain)
    scp.logg.info(code)

    return None
