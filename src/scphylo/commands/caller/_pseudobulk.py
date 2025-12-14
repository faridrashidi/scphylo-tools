import glob
import subprocess

import click
import pandas as pd

import scphylo as scp
from scphylo.ul._servers import cmd, write_cmds_get_main


@click.command(short_help="Run Pseudo Bulk Calling.")
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
    default=["2:00:00", "10:00:00"],
    multiple=True,
    show_default=True,
    help="Time List [indexing, mapping].",
)
@click.option(
    "--mem",
    default=["50", "50"],
    multiple=True,
    show_default=True,
    help="Memory List [indexing, mapping].",
)
@click.option(
    "--afterok",
    default=None,
    type=str,
    show_default=True,
    help="Afterok.",
)
def pseudobulk(outdir, ref, time, mem, afterok):
    """Run Pseudo Bulk Calling.

    scphylo caller pseudobulk ...

    BAM files --> VCF file
    """
    if ref == "hg19":
        config = scp.settings.hg19
    elif ref == "mm10":
        config = scp.settings.mm10
    elif ref == "hg38":
        config = scp.settings.hg38
    else:
        config = None

    if ref == "hg19" or ref == "hg38":
        chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
    elif ref == "mm10":
        chroms = [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY", "chrM"]

    def step1(afterok):
        def get_command(sample):
            cmds = ""
            cmds += cmd([f"mkdir -p {outdir}/_indexing/pseudobulk"])
            cmds += cmd([f"module load {scp.settings.tools['samtools']}"])
            cmds += cmd(
                [
                    "samtools merge",
                    "-in",
                    f"{' -in '.join([])}",
                    f"-out {outdir}/_indexing/pseudobulk/merged.bam",
                ]
            )
            cmds += cmd(["echo Done!"], islast=True)
            return cmds

        df_cmds["cmd"] = df_cmds.apply(lambda x: get_command(x["sample"]), axis=1)
        cmdmain = write_cmds_get_main(
            df_cmds,
            "hapcaller-1of3",
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
        def get_command(chrom):
            files = " --variant ".join(
                [f"{outdir}/{s}.hapcaller.g.vcf" for s in df_cmds["sample"].values]
            )
            cmds = ""
            cmds += cmd([f"module load {scp.settings.tools['gatk']}"])
            cmds += cmd([f"mkdir -p {outdir}/_calling"])
            cmds += cmd(
                [
                    f'gatk --java-options "-Xmx{int(mem[1])-10}g"',
                    "GenotypeGVCFs",
                    f"--reference {config['ref']}",
                    f"--variant {files}",
                    f"--output {outdir}/_calling/jointcalls.{chrom}.g.vcf",
                    # f"-nt 32",
                    f"--intervals {chrom}",
                ]
            )
            cmds += cmd(["echo Done!"], islast=True)
            return cmds

        df_cmds_tmp = pd.DataFrame()
        df_cmds_tmp["chrom"] = chroms
        df_cmds_tmp["cmd"] = df_cmds_tmp.apply(
            lambda x: get_command(x["chrom"]), axis=1
        )
        cmdmain = write_cmds_get_main(
            df_cmds_tmp,
            "hapcaller-2of3",
            time[1],
            mem[1],
            None,
            1,
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

    cmd1 = step1(afterok)
    code1 = subprocess.getoutput(cmd1)
    cmd2 = step2(code1)
    code2 = subprocess.getoutput(cmd2)
    scp.logg.info(code2)

    return None
