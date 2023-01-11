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

    scphylo caller varscan2 /path/to/in/dir normal_name hg19|hg38|mm10

    BAM files (*.markdup_bqsr.bam) --> VCF files (*.varscan2.vcf)
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
            "varscan2-1of2",
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
            cmds += cmd([f"mkdir -p {outdir}/{sample}.varscan2"])
            cmds += cmd(
                [
                    "varscan",
                    "somatic",
                    f"{outdir}/{normal}.markdup_bqsr.pileup",
                    f"{outdir}/{sample}.markdup_bqsr.pileup",
                    f"--output-snp {outdir}/{sample}.varscan2/snp.vcf",
                    f"--output-indel {outdir}/{sample}.varscan2/indel.vcf",
                    "--output-vcf",
                ]
            )
            cmds += cmd([f"module load {scp.settings.tools['bcftools']}"])
            cmds += cmd(
                [
                    "bgzip",
                    f"{outdir}/{sample}.varscan2/snp.vcf",
                ]
            )
            cmds += cmd(
                [
                    "bgzip",
                    f"{outdir}/{sample}.varscan2/indel.vcf",
                ]
            )
            cmds += cmd(
                [
                    "bcftools",
                    "index",
                    f"{outdir}/{sample}.varscan2/snp.vcf.gz",
                ]
            )
            cmds += cmd(
                [
                    "bcftools",
                    "index",
                    f"{outdir}/{sample}.varscan2/indel.vcf.gz",
                ]
            )
            cmds += cmd(
                [
                    "bcftools",
                    "concat",
                    f"{outdir}/{sample}.varscan2/snp.vcf.gz",
                    f"{outdir}/{sample}.varscan2/indel.vcf.gz",
                    f"--output {outdir}/{sample}.varscan2.tmp.vcf",
                    "--output-type v",
                    "--allow-overlaps",
                ]
            )

            cmds += cmd(
                [
                    "varscan",
                    "processSomatic",
                    f"{outdir}/{sample}.varscan2.tmp.vcf",
                ]
            )
            cmds += cmd(
                [
                    "bgzip",
                    f"{outdir}/{sample}.varscan2.tmp.Somatic.vcf",
                ]
            )
            cmds += cmd(
                [
                    "bgzip",
                    f"{outdir}/{sample}.varscan2.tmp.LOH.vcf",
                ]
            )
            cmds += cmd(
                [
                    "bcftools",
                    "index",
                    f"{outdir}/{sample}.varscan2.tmp.Somatic.vcf.gz",
                ]
            )
            cmds += cmd(
                [
                    "bcftools",
                    "index",
                    f"{outdir}/{sample}.varscan2.tmp.LOH.vcf.gz",
                ]
            )
            cmds += cmd(
                [
                    "bcftools",
                    "concat",
                    f"{outdir}/{sample}.varscan2.tmp.Somatic.vcf.gz",
                    f"{outdir}/{sample}.varscan2.tmp.LOH.vcf.gz",
                    f"--output {outdir}/{sample}.varscan2.vcf",
                    "--output-type v",
                    "--allow-overlaps",
                ]
            )
            cmds += cmd(
                [
                    "bcftools",
                    "sort",
                    f"{outdir}/{sample}.varscan2.vcf",
                    f"--output-file {outdir}/{sample}.varscan2.vcf",
                    "--output-type v",
                ]
            )

            cmds += cmd(
                [
                    f"rm -rf {outdir}/{sample}.varscan2",
                    f"rm -rf {outdir}/{sample}.markdup_bqsr.pileup",
                    f"rm -rf {outdir}/{sample}.varscan2.tmp.Germline.hc.vcf",
                    f"rm -rf {outdir}/{sample}.varscan2.tmp.Germline.vcf",
                    f"rm -rf {outdir}/{sample}.varscan2.tmp.Somatic.hc.vcf",
                    f"rm -rf {outdir}/{sample}.varscan2.tmp.LOH.hc.vcf",
                    f"rm -rf {outdir}/{sample}.varscan2.tmp.vcf",
                    f"rm -rf {outdir}/{sample}.varscan2.tmp.Somatic.vcf.gz",
                    f"rm -rf {outdir}/{sample}.varscan2.tmp.Somatic.vcf.gz.csi",
                    f"rm -rf {outdir}/{sample}.varscan2.tmp.LOH.vcf.gz",
                    f"rm -rf {outdir}/{sample}.varscan2.tmp.LOH.vcf.gz.csi",
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
            "varscan2-2of2",
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
    scp.logg.info(code2)

    return None
