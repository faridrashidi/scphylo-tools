import os

import click
import pandas as pd

import scphylo as scp
from scphylo.ul._servers import cmd, write_cmds_get_main


@click.command(short_help="Run vartrix.")
@click.argument(
    "bam",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument(
    "barcodes",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument(
    "vcf",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument(
    "ref",
    required=True,
    type=click.Choice(scp.settings.refs),
)
@click.option(
    "--time",
    default="2-00:00:00",
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
def vartrix(bam, barcodes, vcf, ref, time, mem, afterok):
    """Run vartrix.

    scphylo caller vartrix /path/bam /path/barcodes /path/vcf hg19|hg38|mm10|m10x

    * --> MTX matrices
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

    outdir, _ = scp.ul.dir_base(bam)
    scp.ul.mkdir(f"{outdir}/_vartrix")

    def cmds():
        cmds = ""
        cmds += cmd([f"module load {scp.settings.tools['vartrix']}"])
        cmds += cmd([f"cp {barcodes} {outdir}/_vartrix/"])
        cmds += cmd(
            [
                "vartrix",
                f"--bam {bam}",
                f"--cell-barcodes {barcodes}",
                f"--vcf {vcf}",
                f"--fasta {config['ref']}",
                "--threads $SLURM_CPUS_PER_TASK",
                "--scoring-method coverage",  # coverage, consensus, alt_frac
                "--umi",
                "--no-duplicates",
                f"--out-matrix {outdir}/_vartrix/alt.mtx",
                f"--ref-matrix {outdir}/_vartrix/ref.mtx",
                f"--out-variants {outdir}/_vartrix/var_names.txt"
                # f"--out-matrix {outdir}/_vartrix/gen.mtx",
            ]
        )
        cmds += cmd(["echo Done!"], islast=True)
        return cmds

    df_cmds = pd.DataFrame()
    df_cmds["cmd"] = [cmds()]

    cmdmain = write_cmds_get_main(
        df_cmds,
        "vartrix",
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
