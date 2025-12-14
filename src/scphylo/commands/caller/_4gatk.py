import glob
import subprocess

import click
import pandas as pd

import scphylo as scp
from scphylo.ul._servers import cmd, write_cmds_get_main


@click.command(short_help="Run GATK.")
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
@click.argument(
    "dna_or_rna",
    required=True,
    type=click.Choice(["dna", "rna"]),
)
@click.option(
    "--time",
    default="0-20:00:00",
    type=str,
    show_default=True,
    help="Time.",
)
@click.option(
    "--mem",
    default="90",
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
def gatk(outdir, ref, dna_or_rna, time, mem, afterok):
    """Run GATK.

    caller caller gatk /path/to/in/dir hg19|hg38|mm10 dna|rna

    BAM files (*.mapped.bam) --> BAM files (*.markdup_bqsr.bam)
    """
    if ref == "hg19":
        config = scp.settings.hg19
    elif ref == "mm10":
        config = scp.settings.mm10
    elif ref == "hg38":
        config = scp.settings.hg38
    else:
        config = None

    def get_command(sample):
        cmds = ""
        cmds += cmd([f"module load {scp.settings.tools['gatk']}"])
        cmds += cmd(
            [
                f'gatk --java-options "-Xmx{int(mem)-10}g"',
                "AddOrReplaceReadGroups",
                f"--INPUT {outdir}/{sample}.mapped.bam",
                f"--OUTPUT {outdir}/{sample}.rgadded.bam",
                "--SORT_ORDER coordinate",
                f"--RGID {sample}",
                "--RGLB lib1",
                "--RGPL ILLUMINA",
                "--RGPU unit1",
                f"--RGSM {sample}",
            ]
        )
        cmds += cmd(
            [
                f'gatk --java-options "-Xmx{int(mem)-10}g"',
                "MarkDuplicates",
                f"--INPUT {outdir}/{sample}.rgadded.bam",
                f"--OUTPUT {outdir}/{sample}.markdup.bam",
                f"--METRICS_FILE {outdir}/{sample}.metrics",
            ]
        )
        cmds += cmd(
            [
                "rm -rf",
                f"{outdir}/{sample}.rgadded.bam",
                f"{outdir}/{sample}.metrics",
            ]
        )
        if dna_or_rna == "rna":
            cmds += cmd(
                [
                    f'gatk --java-options "-Xmx{int(mem)-10}g"',
                    "SplitNCigarReads",
                    f"--reference {config['ref']}",
                    f"--input {outdir}/{sample}.markdup.bam",
                    f"--output {outdir}/{sample}.splitted.bam",
                ]
            )
        cmds += cmd(
            [
                f'gatk --java-options "-Xmx{int(mem)-10}g"',
                "BaseRecalibrator",
                f"--reference {config['ref']}",
                f"--input {outdir}/{sample}.splitted.bam"
                if dna_or_rna == "rna"
                else f"--input {outdir}/{sample}.markdup.bam",
                f"--output {outdir}/{sample}.markdup_bqsr.report",
                "--known-sites",
                f"{' --known-sites '.join(config['known_sites'])}",
            ]
        )
        cmds += cmd(
            [
                f'gatk --java-options "-Xmx{int(mem)-10}g"',
                "ApplyBQSR",
                f"--reference {config['ref']}",
                f"--bqsr-recal-file {outdir}/{sample}.markdup_bqsr.report",
                f"--input {outdir}/{sample}.splitted.bam"
                if dna_or_rna == "rna"
                else f"--input {outdir}/{sample}.markdup.bam",
                f"--output {outdir}/{sample}.markdup_bqsr.bam",
            ]
        )
        cmds += cmd(
            [
                "rm -rf",
                f"{outdir}/{sample}.markdup.bam",
                f"{outdir}/{sample}.markdup_bqsr.report",
            ]
        )
        if dna_or_rna == "rna":
            cmds += cmd(
                [
                    "rm -rf",
                    f"{outdir}/{sample}.splitted.bam",
                    f"{outdir}/{sample}.splitted.bai",
                ]
            )
        cmds += cmd(["echo Done!"], islast=True)
        return cmds

    files = glob.glob(f"{outdir}/*.mapped.bam")
    temp = []
    for file in files:
        file = file.split("/")[-1].replace(".mapped.bam", "")
        temp.append({"sample": file})
    df_cmds = pd.DataFrame(temp)
    df_cmds["cmd"] = df_cmds.apply(lambda x: get_command(x["sample"]), axis=1)

    cmdmain = write_cmds_get_main(
        df_cmds,
        "gatk",
        time,
        mem,
        None,
        1,
        scp.settings.tools["email"],
        f"{outdir}/_tmp",
        afterok,
    )
    code = subprocess.getoutput(cmdmain)
    scp.logg.info(code)

    return None
