import subprocess

import click

import scphylo as scp
from scphylo.ul._servers import cmd, write_cmds_get_main


@click.command(short_help="Run VarScan2.")
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
def varscan2(outdir, normal, ref, time, mem, afterok):
    """Run VarScan2.

    scphylo varscan2 /path/to/in/dir normal_name hg19|hg38|mm10

    BAM files (*.markdup_bqsr.bam) --> VCF files (*.varscan2.vcf)
    """
    if ref == "hg19":
        config = scp.settings.hg19
    elif ref == "mm10":
        config = scp.settings.mm10
    elif ref == "hg38":
        config = scp.settings.hg38

    def step1(afterok):
        def get_command(sample):
            cmds = ""
            cmds += cmd([f"module load {scp.settings.tools['samtools']}"])
            cmds += cmd(
                [
                    "samtools",
                    "mpileup",
                    "-B",
                    f"-f {config['ref']}",
                    f"{outdir}/{sample}.markdup_bqsr.bam",
                    f"> {outdir}/{sample}.markdup_bqsr.pileup",
                ]
            )
            cmds += cmd(["echo Done!"], islast=True)
            return cmds

        df_cmds_tmp = scp.ul.get_samples_df(outdir, normal=None)
        df_cmds_tmp["cmd"] = df_cmds_tmp.apply(
            lambda x: get_command(x["sample"]), axis=1
        )
        cmdmain = write_cmds_get_main(
            df_cmds_tmp,
            "varscan2-1of3",
            time,
            mem,
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
            cmds += cmd([f"module load {scp.settings.tools['varscan']}"])
            cmds += cmd(
                [
                    "varscan",
                    "somatic",
                    f"{outdir}/{normal}.markdup_bqsr.pileup",
                    f"{outdir}/{sample}.markdup_bqsr.pileup",
                    "--output-vcf",
                ]
            )
            cmds += cmd(["echo Done!"], islast=True)
            return cmds

        df_cmds_tmp = scp.ul.get_samples_df(outdir, normal)
        df_cmds_tmp["cmd"] = df_cmds_tmp.apply(
            lambda x: get_command(x["sample"]), axis=1
        )
        cmdmain = write_cmds_get_main(
            df_cmds_tmp,
            "varscan2-2of3",
            time,
            mem,
            None,
            1,
            scp.settings.tools["email"],
            f"{outdir}/_tmp",
            afterok,
        )
        return cmdmain

    def step3(afterok):
        def get_command(sample):
            cmds = ""
            cmds += cmd(
                [
                    f"rm -rf {outdir}/{sample}.markdup_bqsr.pileup",
                ]
            )
            cmds += cmd(["echo Done!"], islast=True)
            return cmds

        df_cmds_tmp = scp.ul.get_samples_df(outdir, normal=None)
        df_cmds_tmp["cmd"] = df_cmds_tmp.apply(
            lambda x: get_command(x["sample"]), axis=1
        )
        cmdmain = write_cmds_get_main(
            df_cmds_tmp,
            "varscan2-3of3",
            time,
            mem,
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
    cmd3 = step3(code2)
    code3 = subprocess.getoutput(cmd3)
    scp.logg.info(code3)

    return None
