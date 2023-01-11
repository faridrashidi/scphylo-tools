import glob
import subprocess

import click
import pandas as pd

import scphylo as scp
from scphylo.ul._servers import cmd, write_cmds_get_main


@click.command(short_help="Run bam-readcount.")
@click.argument(
    "outdir",
    required=True,
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, readable=True, resolve_path=True
    ),
)
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
    default=["20:00:00", "01:00:00"],
    multiple=True,
    show_default=True,
    help="Time List [mapping, collecting].",
)
@click.option(
    "--mem",
    default=["20", "20"],
    multiple=True,
    show_default=True,
    help="Memory List [mapping, collecting].",
)
@click.option(
    "--afterok",
    default=None,
    type=str,
    show_default=True,
    help="Afterok.",
)
def bamreadcount(outdir, infile, ref, time, mem, afterok):
    """Run bamreadcount.

    scphylo caller bamreadcount /path/to/bam/dir /path/to/snv/file hg19|hg38|mm10

    The snv file must contains `CHROM,POS,Allele`

    BAM files (*.markdup_bqsr.bam) --> AnnData file (_bamreadcount.h5ad.gz)
    """
    if ref == "hg19":
        config = scp.settings.hg19
    elif ref == "mm10":
        config = scp.settings.mm10
    elif ref == "hg38":
        config = scp.settings.hg38
    else:
        config = None

    def step1(afterok):
        def get_command(sample):
            cmds = ""
            cmds += cmd([f"module load {scp.settings.tools['bamreadcount']}"])
            cmds += cmd(
                [
                    "bam-readcount",
                    "-w 0",
                    f"-f {config['ref']}",
                    f"-l {outdir}/_tmp/bamreadcount.input",
                    f"{outdir}/{sample}.markdup_bqsr.bam",
                    f"> {outdir}/{sample}.bamreadcount",
                ]
            )
            cmds += cmd(["echo Done!"], islast=True)
            return cmds

        df_cmds["cmd"] = df_cmds.apply(lambda x: get_command(x["sample"]), axis=1)
        cmdmain = write_cmds_get_main(
            df_cmds,
            "bamreadcount-1of2",
            time[0],
            mem[0],
            None,
            1,
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
                    "python "
                    f"{scp.ul.get_file('scphylo.ul/bamreadcount.py')} "
                    f"{outdir} "
                    f"{infile}"
                ]
            )
            cmds += cmd(["echo Done!"], islast=True)
            return cmds

        cmds_tmp = pd.DataFrame()
        cmds_tmp["cmd"] = [get_command()]
        cmdmain = write_cmds_get_main(
            cmds_tmp,
            "bamreadcount-2of2",
            time[1],
            mem[1],
            None,
            2,
            scp.settings.tools["email"],
            f"{outdir}/_tmp",
            afterok,
        )
        return cmdmain

    files = glob.glob(f"{outdir}/*.markdup_bqsr.bam")
    temp = []
    for file in files:
        file = file.split("/")[-1].replace(".markdup_bqsr.bam", "")
        temp.append({"sample": file})
    df_cmds = pd.DataFrame(temp)

    scp.ul.mkdir(f"{outdir}/_tmp")
    vcf = pd.read_csv(infile)
    vcf[["CHROM", "POS", "POS"]].to_csv(
        f"{outdir}/_tmp/bamreadcount.input", sep="\t", header=None, index=None
    )

    cmd1 = step1(afterok)
    code1 = subprocess.getoutput(cmd1)
    cmd2 = step2(code1)
    code2 = subprocess.getoutput(cmd2)
    scp.logg.info(code2)

    return None
