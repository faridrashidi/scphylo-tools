import subprocess

import click
import pandas as pd

import scphylo as scp
from scphylo.ul._servers import cmd, write_cmds_get_main


@click.command(short_help="Run SnpEff.")
@click.argument(
    "infile",
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
def snpeff(infile, ref, time, mem, afterok):
    """Run SnpEff.

    scphylo snpeff /path/to/vcf/file hg19|hg38|mm10

    VCF file --> AnnData file (_genotype.h5ad.gz)
    """
    outdir, _ = scp.ul.dir_base(infile)

    def step1(afterok):
        def get_command():
            cmds = ""
            cmds += cmd([f"module load {scp.settings.tools['snpeff']}"])
            cmds += cmd(
                [
                    f"java -Xmx{int(mem)-4}g -jar $SNPEFF_JAR",
                    "-hgvs",
                    # "-cancer",
                    # f"-cancerSamples {1}",
                    ref,
                    infile,
                    f"-stats {infile[:-len('.vcf')] + '.snpeff'}"
                    f"> {infile[:-len('.vcf')] + '.ann.vcf'}",
                ]
            )
            cmds += cmd(["rm -rf", f"{infile[:-len('.vcf')] + '.snpeff.genes.txt'}"])
            cmds += cmd(
                [
                    "mv",
                    f"{infile[:-len('.vcf')] + '.snpeff'}",
                    f"{infile[:-len('.vcf')] + '.snpeff.html'}",
                ]
            )
            cmds += cmd(["echo Done!"], islast=True)
            return cmds

        df_cmds = pd.DataFrame()
        df_cmds["cmd"] = [get_command()]
        cmdmain = write_cmds_get_main(
            df_cmds,
            "snpeff-1of2",
            time,
            mem,
            None,
            2,
            scp.settings.tools["email"],
            f"{outdir}/_tmp",
            afterok,
        )
        return cmdmain

    def step2(afterok):
        def get_command():
            cmds = ""
            cmds += cmd(
                [
                    f"python {scp.ul.get_file('scphylo.ul/cyvcf.py')}",
                    f"{infile[:-len('.vcf')] + '.ann.vcf'}",
                ]
            )
            cmds += cmd(["echo Done!"], islast=True)
            return cmds

        cmds_tmp = pd.DataFrame()
        cmds_tmp["cmd"] = [get_command()]
        cmdmain = write_cmds_get_main(
            cmds_tmp,
            "snpeff-2of2",
            time,
            mem,
            None,
            2,
            scp.settings.tools["email"],
            f"{outdir}/_tmp",
            afterok,
        )
        return cmdmain

    cmd1 = step1(afterok)
    code1 = subprocess.getoutput(cmd1)
    cmd2 = step2(code1)
    code2 = subprocess.getoutput(cmd2)
    scp.logg.info(code2)

    return None
