import subprocess

import click

import scphylo as scp
from scphylo.ul._servers import cmd, write_cmds_get_main


@click.command(short_help="Run MuTect2.")
@click.argument(
    "outdir",
    required=True,
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, readable=True, resolve_path=True
    ),
)
@click.argument(
    "normal",
    required=True,
    type=str,
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
    default="50",
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
def mutect2(outdir, normal, ref, time, mem, afterok):
    """Run MuTect2.

    scphylo caller mutect2 /path/to/in/dir normal_name hg19|hg38|mm10

    BAM files (*.markdup_bqsr.bam) --> VCF files (*.mutect2.vcf)
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
                "Mutect2",
                f"--reference {config['ref']}",
                f"--input {outdir}/{sample}.markdup_bqsr.bam",
                f"--input {outdir}/{normal}.markdup_bqsr.bam",
                f"--normal-sample {normal}",
                f"--output {outdir}/{sample}.tmp.mutect2.vcf",
            ]
        )
        cmds += cmd(
            [
                f'gatk --java-options "-Xmx{int(mem)-10}g"',
                "FilterMutectCalls",
                f"--reference {config['ref']}",
                f"--variant {outdir}/{sample}.tmp.mutect2.vcf",
                f"--output {outdir}/{sample}.mutect2.vcf",
            ]
        )
        cmds += cmd(
            [
                "rm -rf",
                f"{outdir}/{sample}.mutect2.vcf.filteringStats.tsv",
                f"{outdir}/{sample}.mutect2.vcf.stats",
                f"{outdir}/{sample}.mutect2.vcf.idx",
                f"{outdir}/{sample}.tmp.mutect2.vcf",
                f"{outdir}/{sample}.tmp.mutect2.vcf.stats",
                f"{outdir}/{sample}.tmp.mutect2.vcf.idx",
            ]
        )
        cmds += cmd(["echo Done!"], islast=True)
        return cmds

    df_cmds = scp.ul.get_samples_df(outdir, normal)
    df_cmds["cmd"] = df_cmds.apply(lambda x: get_command(x["sample"]), axis=1)

    cmdmain = write_cmds_get_main(
        df_cmds,
        "mutect2",
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
