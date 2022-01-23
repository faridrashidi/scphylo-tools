import subprocess

import click

import scphylo as scp
from scphylo.ul._servers import cmd, write_cmds_get_main


@click.command(short_help="Run Strelka.")
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
def strelka(outdir, normal, ref, time, mem, afterok):
    """Run Strelka.

    scphylo strelka /path/to/in/dir normal_name hg19|hg38|mm10

    BAM files (*.markdup_bqsr.bam) --> VCF files (*.strelka.vcf)
    """
    if ref == "hg19":
        config = scp.settings.hg19
    elif ref == "mm10":
        config = scp.settings.mm10
    elif ref == "hg38":
        config = scp.settings.hg38

    def get_command(sample):
        cmds = ""
        cmds += cmd([f"module load {scp.settings.tools['strelka']}"])
        cmds += cmd(
            [
                "cofigureStrelkaSomaticWorkflow.py",
                f"--referenceFasta={config['ref']}",
                f"--normalBam={outdir}/{normal}.markdup_bqsr.bam",
                f"--tumorBam={outdir}/{sample}.markdup_bqsr.bam",
                f"--runDir={outdir}",
            ]
        )
        cmds += cmd(["demo_out/runWorkflow.py -m local -j $SLURM_CPUS_PER_TASK"])
        cmds += cmd(["echo Done!"], islast=True)
        return cmds

    df_cmds = scp.ul.get_samples_df(outdir, normal)
    df_cmds["cmd"] = df_cmds.apply(lambda x: get_command(x["sample"]), axis=1)

    cmdmain = write_cmds_get_main(
        df_cmds,
        "strelka",
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
