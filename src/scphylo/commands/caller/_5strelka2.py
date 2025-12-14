import subprocess

import click

import scphylo as scp
from scphylo.ul._servers import cmd, write_cmds_get_main


@click.command(short_help="Run Strelka2.")
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
def strelka2(outdir, normal, ref, time, mem, afterok):
    """Run Strelka2.

    scphylo caller strelka2 /path/to/in/dir normal_name hg19|hg38|mm10

    BAM files (*.markdup_bqsr.bam) --> VCF files (*.strelka2.vcf)
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
        cmds += cmd([f"module load {scp.settings.tools['strelka']}"])
        cmds += cmd(
            [
                "configureStrelkaSomaticWorkflow.py",
                f"--referenceFasta={config['ref']}",
                f"--normalBam={outdir}/{normal}.markdup_bqsr.bam",
                f"--tumorBam={outdir}/{sample}.markdup_bqsr.bam",
                f"--runDir={outdir}/{sample}.strelka2",
            ]
        )
        cmds += cmd(
            [
                f"{outdir}/{sample}.strelka2/runWorkflow.py -m local -j "
                "$SLURM_CPUS_PER_TASK"
            ]
        )
        cmds += cmd([f"module load {scp.settings.tools['bcftools']}"])
        cmds += cmd(
            [
                "bcftools",
                "concat",
                f"{outdir}/{sample}.strelka2/results/variants/somatic.snvs.vcf.gz",
                f"{outdir}/{sample}.strelka2/results/variants/somatic.indels.vcf.gz",
                f"--output {outdir}/{sample}.strelka2.tmp.vcf",
                "--output-type v",
                "--allow-overlaps",
            ]
        )
        cmds += cmd(
            [
                "bcftools",
                "view --apply-filters PASS",
                f"{outdir}/{sample}.strelka2.tmp.vcf",
                f"> {outdir}/{sample}.strelka2.vcf",
            ]
        )

        cmds += cmd(
            [
                "bcftools",
                "sort",
                f"{outdir}/{sample}.strelka2.vcf",
                f"--output-file {outdir}/{sample}.strelka2.vcf",
                "--output-type v",
            ]
        )
        cmds += cmd(
            [
                f"rm -rf {outdir}/{sample}.strelka2",
                f"rm -rf {outdir}/{sample}.strelka2.tmp.vcf",
            ]
        )
        cmds += cmd(["echo Done!"], islast=True)
        return cmds

    df_cmds = scp.ul.get_samples_df(outdir, normal)
    df_cmds["cmd"] = df_cmds.apply(lambda x: get_command(x["sample"]), axis=1)

    cmdmain = write_cmds_get_main(
        df_cmds,
        "strelka2",
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
