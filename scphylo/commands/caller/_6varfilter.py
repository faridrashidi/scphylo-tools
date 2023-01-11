import glob
import subprocess

import click
import pandas as pd

import scphylo as scp
from scphylo.ul._servers import cmd, write_cmds_get_main


@click.command(short_help="Run VariantFiltration.")
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
    default="1-00:00:00",
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
def varfilter(outdir, ref, time, mem, afterok):
    """Run varfilter.

    scphylo caller filter /path/to/vcf/file hg19|hg38|mm10

    VCF file --> VCF.filtered file
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
                "IndexFeatureFile",
                f"--input {outdir}/{sample}.vcf",
            ]
        )
        cmds += cmd(
            [
                f'gatk --java-options "-Xmx{int(mem)-10}g"',
                "VariantFiltration",
                f"--reference {config['ref']}",
                f"--variant {outdir}/{sample}.vcf",
                f"--output {outdir}/{sample}.tmp.filtered.vcf",
                "--cluster-window-size 35",
                "--cluster-size 3",
                '--filter-name "FS"',
                '--filter-expression "FS > 30.0"',
                '--filter-name "QD"',
                '--filter-expression "QD < 2.0"',
            ]
        )
        cmds += cmd(
            [
                f'gatk --java-options "-Xmx{int(mem)-10}g"',
                "SelectVariants",
                f"--variant {outdir}/{sample}.tmp.filtered.vcf",
                f"--output {outdir}/{sample}.filtered.vcf",
                "--exclude-filtered",
            ]
        )
        cmds += cmd(
            [
                "rm -rf",
                f"{outdir}/{sample}.tmp.filtered.vcf",
            ]
        )
        cmds += cmd(["echo Done!"], islast=True)
        return cmds

    files = glob.glob(f"{outdir}/*.vcf")
    temp = []
    for file in files:
        file = file.split("/")[-1].replace(".vcf", "")
        temp.append({"sample": file})
    df_cmds = pd.DataFrame(temp)
    df_cmds["cmd"] = df_cmds.apply(lambda x: get_command(x["sample"]), axis=1)

    cmdmain = write_cmds_get_main(
        df_cmds,
        "varfilter",
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
