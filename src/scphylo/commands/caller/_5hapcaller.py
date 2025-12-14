import glob
import subprocess

import click
import pandas as pd

import scphylo as scp
from scphylo.ul._servers import cmd, write_cmds_get_main


@click.command(short_help="Run HaplotypeCaller.")
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
    "dna_or_rna",
    required=True,
    type=click.Choice(["dna", "rna"]),
)
@click.option(
    "--time",
    default=["20:00:00", "2-00:00:00"],
    multiple=True,
    show_default=True,
    help="Time List [indexing, mapping].",
)
@click.option(
    "--mem",
    default=["90", "90"],
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
def hapcaller(outdir, ref, dna_or_rna, time, mem, afterok):
    """Run HaplotypeCaller.

    scphylo caller hapcaller /path/to/in/dir hg19|hg38|mm10 dna|rna

    BAM files (*.markdup_bqsr.bam) --> VCF file (_hapcaller/jointcalls.g.vcf)
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
            cmds += cmd([f"module load {scp.settings.tools['gatk']}"])
            if dna_or_rna == "rna":
                cmds += cmd(
                    [
                        f'gatk --java-options "-Xmx{int(mem[0])-10}g"',
                        "HaplotypeCaller",
                        f"--reference {config['ref']}",
                        f"--input {outdir}/{sample}.markdup_bqsr.bam",
                        f"--output {outdir}/{sample}.hapcaller.g.vcf",
                        "--dont-use-soft-clipped-bases true",
                        "--standard-min-confidence-threshold-for-calling 20",
                        "--emit-ref-confidence GVCF",
                        "--disable-read-filter MappingQualityReadFilter",
                        "--disable-read-filter GoodCigarReadFilter",
                        "--disable-read-filter NotSecondaryAlignmentReadFilter",
                        "--disable-read-filter MappedReadFilter",
                        "--disable-read-filter MappingQualityAvailableReadFilter",
                        "--disable-read-filter",
                        "NonZeroReferenceLengthAlignmentReadFilter",
                        "--disable-read-filter NotDuplicateReadFilter",
                        "--disable-read-filter PassesVendorQualityCheckReadFilter",
                        "--disable-read-filter WellformedReadFilter",
                    ]
                )
            else:
                cmds += cmd(
                    [
                        f'gatk --java-options "-Xmx{int(mem[0])-10}g"',
                        "HaplotypeCaller",
                        f"--reference {config['ref']}",
                        f"--input {outdir}/{sample}.markdup_bqsr.bam",
                        f"--output {outdir}/{sample}.hapcaller.g.vcf",
                        "--dont-use-soft-clipped-bases true",
                        "--standard-min-confidence-threshold-for-calling 20",
                        "--emit-ref-confidence GVCF",
                    ]
                )
            cmds += cmd(["echo Done!"], islast=True)
            return cmds

        df_cmds["cmd"] = df_cmds.apply(lambda x: get_command(x["sample"]), axis=1)
        cmdmain = write_cmds_get_main(
            df_cmds,
            "hapcaller-1of4",
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
            cmds = ""
            files = " --variant ".join(
                [f"{outdir}/{s}.hapcaller.g.vcf" for s in df_cmds["sample"].values]
            )
            cmds += cmd([f"module load {scp.settings.tools['gatk']}"])
            cmds += cmd([f"mkdir -p {outdir}/_hapcaller"])
            cmds += cmd(
                [
                    f'gatk --java-options "-Xmx{int(mem[1])-10}g"',
                    "CombineGVCFs",
                    f"--reference {config['ref']}",
                    f"--variant {files}",
                    f"--output {outdir}/_hapcaller/combinedcalls.{chrom}.g.vcf",
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
            "hapcaller-2of4",
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
        def get_command(chrom):
            cmds = ""
            cmds += cmd([f"module load {scp.settings.tools['gatk']}"])
            cmds += cmd([f"mkdir -p {outdir}/_hapcaller"])
            cmds += cmd(
                [
                    f'gatk --java-options "-Xmx{int(mem[1])-10}g"',
                    "GenotypeGVCFs",
                    f"--reference {config['ref']}",
                    f"--variant {outdir}/_hapcaller/combinedcalls.{chrom}.g.vcf",
                    f"--output {outdir}/_hapcaller/jointcalls.{chrom}.g.vcf",
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
            "hapcaller-3of4",
            time[1],
            mem[1],
            None,
            1,
            scp.settings.tools["email"],
            f"{outdir}/_tmp",
            afterok,
        )
        return cmdmain

    def step4(afterok):
        def get_command():
            cmds = ""
            cmds += cmd([f"module load {scp.settings.tools['bcftools']}"])
            files = " ".join(
                [f"{outdir}/_hapcaller/jointcalls.{chrom}.g.vcf" for chrom in chroms]
            )
            cmds += cmd(
                [
                    "bcftools concat",
                    f"-o {outdir}/_hapcaller/jointcalls.g.vcf",
                    f"{files}",
                ]
            )
            cmds += cmd(["echo Done!"], islast=True)
            return cmds

        cmds_tmp = pd.DataFrame()
        cmds_tmp["cmd"] = [get_command()]
        cmdmain = write_cmds_get_main(
            cmds_tmp,
            "hapcaller-4of4",
            time[0],
            mem[0],
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
    cmd3 = step3(code2)
    code3 = subprocess.getoutput(cmd3)
    cmd4 = step4(code3)
    code4 = subprocess.getoutput(cmd4)
    scp.logg.info(code4)

    return None
