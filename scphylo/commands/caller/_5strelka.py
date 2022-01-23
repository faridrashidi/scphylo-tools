import glob
import subprocess

import click
import pandas as pd

import scphylo as scp
from scphylo.ul._servers import cmd, write_cmds_get_main


@click.command(short_help="Run Sterlka.")
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
def sterlka(outdir, normal, ref, time, mem, afterok):
    """Run Sterlka.

    scphylo sterlka /path/to/in/dir /path/to/normal/name hg19|hg38|mm10

    BAM files (*.markdup_bqsr.bam) --> VCF files (*.vcf)
    """
    if ref == "hg19":
        config = scp.settings.hg19
    elif ref == "mm10":
        config = scp.settings.mm10
    elif ref == "hg38":
        config = scp.settings.hg38

    def get_command(sample):
        cmds = ""
        cmds += cmd([f"module load {scp.settings.tools['gatk']}"])
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

    files = glob.glob(f"{outdir}/*.markdup_bqsr.bam")
    temp = []
    for file in files:
        file = file.split("/")[-1].replace(".markdup_bqsr.bam", "")
        if file != normal:
            temp.append({"sample": file})
    df_cmds = pd.DataFrame(temp)
    scp.logg.info(
        f"Tumor samples: {','.join(df_cmds.sample.values)} & Normal sample: {normal}"
    )
    df_cmds["cmd"] = df_cmds.apply(lambda x: get_command(x["sample"]), axis=1)

    cmdmain = write_cmds_get_main(
        df_cmds,
        "sterlka",
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
