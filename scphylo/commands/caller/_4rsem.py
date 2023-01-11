import glob
import subprocess

import click
import pandas as pd

import scphylo as scp
from scphylo.ul._servers import cmd, write_cmds_get_main


@click.command(short_help="Run RSEM.")
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
@click.argument(
    "single_or_paired",
    required=True,
    type=click.Choice(["single", "paired"]),
)
@click.option(
    "--time",
    default=["1:00:00", "10:00:00"],
    multiple=True,
    show_default=True,
    help="Time List [indexing, calling].",
)
@click.option(
    "--mem",
    default=["20", "20"],
    multiple=True,
    show_default=True,
    help="Memory List [indexing, calling].",
)
@click.option(
    "--afterok",
    default=None,
    type=str,
    show_default=True,
    help="Afterok.",
)
def rsem(outdir, ref, single_or_paired, time, mem, afterok):
    """Run RSEM.

    scphylo caller rsem /path/to/in/dir hg19|hg38|mm10 paired

    BAM files (*.transcript.bam) --> AnnData file (_expression.h5ad.gz)
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
        def get_command():
            cmds = ""
            cmds += cmd([f"mkdir -p {outdir}/_indexing/rsem"])
            cmds += cmd([f"module load {scp.settings.tools['rsem']}"])
            cmds += cmd(
                [
                    "rsem-prepare-reference",
                    "--num-threads 32",
                    "--gtf",
                    f"{config['annot']}",
                    f"{config['ref']}",
                    f"{outdir}/_indexing/rsem/rsem",
                ]
            )
            cmds += cmd(["echo Done!"], islast=True)
            return cmds

        cmds_tmp = pd.DataFrame()
        cmds_tmp["cmd"] = [get_command()]
        cmdmain = write_cmds_get_main(
            cmds_tmp,
            "rsem-1of3",
            time[0],
            mem[0],
            None,
            32,
            scp.settings.tools["email"],
            f"{outdir}/_tmp",
            afterok,
        )
        return cmdmain

    def step2(afterok):
        def get_command(sample):
            cmds = ""
            cmds += cmd([f"module load {scp.settings.tools['rsem']}"])
            if single_or_paired == "paired":
                cmds += cmd(
                    [
                        "rsem-calculate-expression",
                        "--bam",
                        "--no-bam-output",
                        "--paired-end",
                        "--estimate-rspd",
                        f"{outdir}/{sample}.transcript.bam",
                        f"{outdir}/_indexing/rsem/rsem",
                        f"{outdir}/{sample}",
                    ]
                )
            else:
                cmds += cmd(
                    [
                        "rsem-calculate-expression",
                        "--bam",
                        "--no-bam-output",
                        "--estimate-rspd",
                        f"{outdir}/{sample}.transcript.bam",
                        f"{outdir}/_indexing/rsem/rsem",
                        f"{outdir}/{sample}",
                    ]
                )
            cmds += cmd(
                [
                    "rm -rf",
                    f"{outdir}/{sample}.stat",
                ]
            )
            cmds += cmd(["echo Done!"], islast=True)
            return cmds

        df_cmds["cmd"] = df_cmds.apply(lambda x: get_command(x["sample"]), axis=1)
        cmdmain = write_cmds_get_main(
            df_cmds,
            "rsem-2of3",
            time[1],
            mem[1],
            None,
            1,
            scp.settings.tools["email"],
            f"{outdir}/_tmp",
            afterok,
        )
        return cmdmain

    def step3(afterok):
        def get_command():
            cmds = ""
            cmds += cmd([f"python {scp.ul.get_file('scphylo.ul/rsem.py')} {outdir}"])
            cmds += cmd(["echo Done!"], islast=True)
            return cmds

        cmds_tmp = pd.DataFrame()
        cmds_tmp["cmd"] = [get_command()]
        cmdmain = write_cmds_get_main(
            cmds_tmp,
            "rsem-3of3",
            time[0],
            mem[0],
            None,
            2,
            scp.settings.tools["email"],
            f"{outdir}/_tmp",
            afterok,
        )
        return cmdmain

    files = glob.glob(f"{outdir}/*.transcript.bam")
    temp = []
    for file in files:
        file = file.split("/")[-1].replace(".transcript.bam", "")
        temp.append({"sample": file})
    df_cmds = pd.DataFrame(temp)

    cmd1 = step1(afterok)
    code1 = subprocess.getoutput(cmd1)
    cmd2 = step2(code1)
    code2 = subprocess.getoutput(cmd2)
    cmd3 = step3(code2)
    code3 = subprocess.getoutput(cmd3)
    scp.logg.info(code3)

    return None
