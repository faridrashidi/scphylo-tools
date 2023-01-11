import os

import click

import scphylo as scp
from scphylo.ul._servers import cmd


@click.command(short_help="Run bamsnap.")
# @click.argument(
#     "outdir",
#     required=True,
#     type=click.Path(
#         exists=True, file_okay=False, dir_okay=True, readable=True, resolve_path=True
#     ),
# )
@click.argument(
    "ref",
    required=True,
    type=click.Choice(scp.settings.refs),
)
@click.argument(
    "bam",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument(
    "pos",
    required=True,
    type=str,
)
def bamsnap(ref, bam, pos):
    """Run bamsnap.

    scphylo caller bamsnap ...

    BAM files (*.markdup_bqsr.bam) --> PNG file (_bamsnap.png)
    """
    if ref == "hg19":
        config = scp.settings.hg19
    elif ref == "mm10":
        config = scp.settings.mm10
    elif ref == "hg38":
        config = scp.settings.hg38
    else:
        config = None

    cmds = ""
    # scp.ul.mkdir(f"{outdir}/_bamsnap")
    cmds += cmd(
        [
            "bamsnap",
            f"-ref {config['ref']}",
            f"-bam {bam}",
            f"-pos {pos}",
            f"-out ./{pos}",
            "-silence",
            "-read_color_deletion 'C8C8C8'",
            "-coverage_vaf 0.1",
        ],
        islast=True,
    )
    print(cmds)
    # cmds += cmd(["rm -rf", f"{1}_bamsnap.log"])

    os.system(cmds)

    return None
