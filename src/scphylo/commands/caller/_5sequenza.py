import subprocess

import click
import pandas as pd

import scphylo as scp
from scphylo.ul._servers import cmd, write_cmds_get_main


@click.command(short_help="Run Sequenza.")
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
    default=["1:00:00", "20:00:00"],
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
def sequenza(outdir, normal, ref, time, mem, afterok):
    """Run Sequenza.

    scphylo caller sequenza /path/to/in/dir normal_name hg19|hg38|mm10

    BAM files (*.markdup_bqsr.bam) --> folder (.sequenza)
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
            cmds += cmd([f"mkdir -p {outdir}/_indexing/sequenza"])
            cmds += cmd([f"module load {scp.settings.tools['sequenza']}"])
            cmds += cmd(
                [
                    "sequenza-utils",
                    "gc_wiggle",
                    "-w 50",
                    f"--fasta {config['ref']}",
                    f"-o {outdir}/_indexing/sequenza/sequenza.wig.gz",
                ]
            )
            cmds += cmd(["echo Done!"], islast=True)
            return cmds

        cmds_tmp = pd.DataFrame()
        cmds_tmp["cmd"] = [get_command()]
        cmdmain = write_cmds_get_main(
            cmds_tmp,
            "sequenza-1of2",
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
        def get_command(sample):
            cmds = ""
            cmds += cmd([f"module load {scp.settings.tools['sequenza']}"])
            cmds += cmd([f"module load {scp.settings.tools['samtools']}"])
            cmds += cmd(
                [
                    "sequenza-utils",
                    "bam2seqz",
                    f"--normal {outdir}/{normal}.markdup_bqsr.bam",
                    f"--tumor {outdir}/{sample}.markdup_bqsr.bam",
                    f"--fasta {config['ref']}",
                    f"-gc {outdir}/_indexing/sequenza/sequenza.wig.gz",
                    f"--output {outdir}/{sample}.seqz.gz",
                    # f"--chromosome {' '.join(chroms)}",
                ]
            )
            cmds += cmd(
                [
                    "sequenza-utils",
                    "seqz_binning",
                    f"--seqz {outdir}/{sample}.seqz.gz",
                    "-w 50",
                    f"-o {outdir}/{sample}.50.seqz.gz",
                ]
            )
            cmds += cmd(
                [
                    "Rscript",
                    f"{scp.ul.get_file('scphylo.ul/sequenza.R')}",
                    outdir,
                    sample,
                ]
            )
            # cmds += cmd(
            #     [
            #         "rm -rf",
            #         f"{outdir}/{sample}.50.seqz.gz",
            #         f"{outdir}/{sample}.50.seqz.gz.tbi",
            #         f"{outdir}/{sample}.seqz.gz",
            #         f"{outdir}/{sample}.seqz.gz.tbi",
            #     ]
            # )
            cmds += cmd(["echo Done!"], islast=True)
            return cmds

        df_cmds_tmp = scp.ul.get_samples_df(outdir, normal)
        df_cmds_tmp["cmd"] = df_cmds_tmp.apply(
            lambda x: get_command(x["sample"]), axis=1
        )
        cmdmain = write_cmds_get_main(
            df_cmds_tmp,
            "sequenza-2of2",
            time[1],
            mem[1],
            None,
            1,
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
