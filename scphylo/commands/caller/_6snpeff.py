import glob
import subprocess

import click
import pandas as pd

import scphylo as scp
from scphylo.ul._servers import cmd, write_cmds_get_main


@click.command(short_help="Run SnpEff.")
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
    default="0-02:00:00",
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
def snpeff(outdir, ref, time, mem, afterok):
    """Run SnpEff.

    scphylo caller snpeff /path/to/vcf/file hg19|hg38|mm10

    VCF file --> AnnData file (_genotype.h5ad.gz)
    """

    def get_command(sample):
        cmds = ""
        cmds += cmd([f"module load {scp.settings.tools['snpeff']}"])
        cmds += cmd(
            [
                f"java -Xmx{int(mem)-4}g -jar $SNPEFF_JAR",
                "-hgvs",
                # "-cancer",
                # f"-cancerSamples {1}",
                ref,
                f"{outdir}/{sample}.vcf",
                f"-stats {outdir}/{sample}.snpeff",
                f"> {outdir}/{sample}.ann.vcf",
            ]
        )
        cmds += cmd(["rm -rf", f"{outdir}/{sample}.snpeff.genes.txt"])
        cmds += cmd(
            [
                "mv",
                f"{outdir}/{sample}.snpeff",
                f"{outdir}/{sample}.snpeff.html",
            ]
        )

        cmds += cmd(
            [
                f"python {scp.ul.get_file('scphylo.ul/cyvcf.py')}",
                f"{outdir}/{sample}.ann.vcf",
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
        "snpeff",
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
