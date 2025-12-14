import gzip
import subprocess

import click
import pandas as pd

import scphylo as scp
from scphylo.ul._servers import cmd, write_cmds_get_main


@click.command(short_help="Run STAR.")
@click.argument(
    "indir",
    required=True,
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, readable=True, resolve_path=True
    ),
)
@click.argument(
    "outdir",
    required=True,
    type=click.Path(file_okay=False, dir_okay=True, readable=True, resolve_path=True),
)
@click.argument(
    "ref",
    required=True,
    type=click.Choice(scp.settings.refs),
)
@click.option(
    "--time",
    default=["2:00:00", "20:00:00"],
    multiple=True,
    show_default=True,
    help="Time List [indexing, mapping].",
)
@click.option(
    "--mem",
    default=["50", "90"],
    multiple=True,
    show_default=True,
    help="Memory List [indexing, mapping].",
)
@click.option(
    "--max_multimapping",
    default="10",
    type=str,
    show_default=True,
    help="`outFilterMultimapNmax` parameter in STAR.",
)
@click.option(
    "--afterok",
    default=None,
    type=str,
    show_default=True,
    help="Afterok.",
)
@click.option(
    "--is_pdx",
    default=False,
    type=bool,
    show_default=True,
    is_flag=True,
    help="Is the the PDX.",
)
def star(indir, outdir, ref, time, mem, max_multimapping, afterok, is_pdx):
    """Run STAR.

    scphylo caller star /path/to/in/dir /path/to/out/dir hg19|hg38|mm10 --is_pdx

    FastQ files (*.fastq.gz) --> BAM files (*.mapped.bam)

    if --is_pdx: FastQ files (*.fastq.gz) --> FastQ files (*.fastq.gz)
    """
    scp.ul.mkdir(outdir)
    if ref == "hg19":
        config = scp.settings.hg19
    elif ref == "mm10":
        config = scp.settings.mm10
    elif ref == "hg38":
        config = scp.settings.hg38
    else:
        config = None

    def get_readlength(filename):
        readlength = None
        with gzip.open(filename, "rb") as fin:
            i = 0
            for line in fin:
                readlength = len(line.strip()) - 1
                if i == 1:
                    break
                i += 1
        return readlength

    def step1(afterok):
        def get_command():
            cmds = ""
            cmds += cmd([f"mkdir -p {outdir}/_indexing/star1"])
            cmds += cmd([f"module load {scp.settings.tools['star']}"])
            cmds += cmd(
                [
                    "STAR",
                    "--runMode genomeGenerate",
                    f"--genomeDir {outdir}/_indexing/star1",
                    f"--genomeFastaFiles {config['ref']}",
                    f"--sjdbGTFfile {config['annot']}",
                    f"--sjdbOverhang {readlength}",
                    "--runThreadN 32",
                ]
            )
            cmds += cmd(["echo Done!"], islast=True)
            return cmds

        cmds_tmp = pd.DataFrame()
        cmds_tmp["cmd"] = [get_command()]
        cmdmain = write_cmds_get_main(
            cmds_tmp,
            "star-1of5",
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
        def get_command(sample, is_paired):
            if is_paired:
                infqs = [f"{indir}/{sample}_{i}.fastq.gz" for i in range(1, 3)]
            else:
                infqs = [f"{indir}/{sample}.fastq.gz"]
            cmds = ""
            cmds += cmd([f"module load {scp.settings.tools['star']}"])
            cmds += cmd(
                [
                    "STAR",
                    "--runMode alignReads",
                    f"--genomeDir {outdir}/_indexing/star1",
                    f"--sjdbGTFfile {config['annot']}",
                    f"--outFilterMultimapNmax {max_multimapping}",
                    "--outSAMunmapped None",
                    "--quantMode TranscriptomeSAM GeneCounts",
                    "--runThreadN 1",
                    f"--sjdbOverhang {readlength}",
                    "--readFilesCommand zcat",
                    f"--outFileNamePrefix {outdir}/{sample}_",
                    "--readFilesIn",
                    f"{' '.join(infqs)}",
                ]
            )
            cmds += cmd(
                [
                    "rm -rf",
                    f"{outdir}/{sample}_Aligned.out.sam",
                ]
            )
            cmds += cmd(["echo Done!"], islast=True)
            return cmds

        df_cmds["cmd"] = df_cmds.apply(
            lambda x: get_command(x["sample"], x["is_paired"]), axis=1
        )
        cmdmain = write_cmds_get_main(
            df_cmds,
            "star-2of5",
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
            files = " ".join(
                [f"{outdir}/{s}_SJ.out.tab" for s in df_cmds["sample"].values]
            )
            cmds = ""
            cmds += cmd([f"mkdir -p {outdir}/_indexing/star2"])
            cmds += cmd([f"module load {scp.settings.tools['star']}"])
            cmds += cmd(
                [
                    "STAR",
                    "--runMode genomeGenerate",
                    f"--genomeDir {outdir}/_indexing/star2",
                    f"--genomeFastaFiles {config['ref']}",
                    f"--sjdbGTFfile {config['annot']}",
                    f"--sjdbFileChrStartEnd {files}",
                    f"--sjdbOverhang {readlength}",
                    "--runThreadN 32",
                ]
            )
            cmds += cmd(["echo Done!"], islast=True)
            return cmds

        cmds_tmp = pd.DataFrame()
        cmds_tmp["cmd"] = [get_command()]
        cmdmain = write_cmds_get_main(
            cmds_tmp,
            "star-3of5",
            time[0],
            mem[0],
            None,
            32,
            scp.settings.tools["email"],
            f"{outdir}/_tmp",
            afterok,
        )
        return cmdmain

    def step4(afterok):
        def get_command(sample, is_paired):
            if is_paired:
                infqs = [f"{indir}/{sample}_{i}.fastq.gz" for i in range(1, 3)]
            else:
                infqs = [f"{indir}/{sample}.fastq.gz"]
            cmds = ""
            cmds += cmd([f"module load {scp.settings.tools['star']}"])
            tmp = [
                "STAR",
                "--runMode alignReads",
                f"--genomeDir {outdir}/_indexing/star2",
                "--readFilesCommand zcat",
                "--readFilesIn",
                f"{' '.join(infqs)}",
                f"--outFileNamePrefix {outdir}/{sample}_",
                "--limitBAMsortRAM 30000000000",
                "--outSAMtype BAM SortedByCoordinate",
                f"--sjdbGTFfile {config['annot']}",
                f"--outFilterMultimapNmax {max_multimapping}",
                "--outSAMunmapped None",
                "--quantMode TranscriptomeSAM GeneCounts",
                "--runThreadN 1",
                f"--sjdbOverhang {readlength}",
            ]
            if is_pdx:
                tmp += ["--outReadsUnmapped Fastx"]
            cmds += cmd(tmp)
            cmds += cmd(
                [
                    "rm -rf",
                    f"{outdir}/{sample}_Aligned.out.sam",
                    f"{outdir}/{sample}_Log.progress.out",
                    f"{outdir}/{sample}_Log.out",
                    f"{outdir}/{sample}_ReadsPerGene.out.tab",
                    f"{outdir}/{sample}_SJ.out.tab",
                    f"{outdir}/{sample}__STARgenome",
                ]
            )
            cmds += cmd(
                [
                    "mv",
                    f"{outdir}/{sample}_Aligned.sortedByCoord.out.bam",
                    f"{outdir}/{sample}.mapped.bam",
                ]
            )
            cmds += cmd(
                [
                    "mv",
                    f"{outdir}/{sample}_Aligned.toTranscriptome.out.bam",
                    f"{outdir}/{sample}.transcript.bam",
                ]
            )
            cmds += cmd(
                [
                    "mv",
                    f"{outdir}/{sample}_Log.final.out",
                    f"{outdir}/{sample}.star.log",
                ]
            )
            if is_pdx:
                cmds += cmd(
                    [
                        "mv",
                        f"{outdir}/{sample}_Unmapped.out.mate1",
                        f"{outdir}/{sample}_1.fastq",
                    ]
                )
                cmds += cmd(
                    [
                        "mv",
                        f"{outdir}/{sample}_Unmapped.out.mate2",
                        f"{outdir}/{sample}_2.fastq",
                    ]
                )
                cmds += cmd(
                    [
                        "gzip",
                        f"{outdir}/{sample}_1.fastq",
                    ]
                )
                cmds += cmd(
                    [
                        "gzip",
                        f"{outdir}/{sample}_2.fastq",
                    ]
                )
                cmds += cmd(
                    [
                        "rm -rf",
                        f"{outdir}/{sample}.mapped.bam",
                        f"{outdir}/{sample}.transcript.bam",
                    ]
                )
            cmds += cmd(["echo Done!"], islast=True)
            return cmds

        df_cmds["cmd"] = df_cmds.apply(
            lambda x: get_command(x["sample"], x["is_paired"]), axis=1
        )
        cmdmain = write_cmds_get_main(
            df_cmds,
            "star-4of5",
            time[1],
            mem[1],
            None,
            1,
            scp.settings.tools["email"],
            f"{outdir}/_tmp",
            afterok,
        )
        return cmdmain

    def step5(afterok):
        def get_command():
            cmds = ""
            cmds += cmd(
                [
                    f"python {scp.ul.get_file('scphylo.ul/star.py')}",
                    f"{outdir}",
                ]
            )
            cmds += cmd(
                [
                    "rm -rf",
                    f"{outdir}/*.star.log",
                ]
            )
            cmds += cmd(["echo Done!"], islast=True)
            return cmds

        cmds_tmp = pd.DataFrame()
        cmds_tmp["cmd"] = [get_command()]
        cmdmain = write_cmds_get_main(
            cmds_tmp,
            "star-5of5",
            time[0],
            mem[0],
            None,
            2,
            scp.settings.tools["email"],
            f"{outdir}/_tmp",
            afterok,
        )
        return cmdmain

    temp = scp.ul.is_paired(indir)
    if temp[0]["is_paired"]:
        readlength = get_readlength(f"{indir}/{temp[0]['sample']}_1.fastq.gz")
    else:
        readlength = get_readlength(f"{indir}/{temp[0]['sample']}.fastq.gz")
    df_cmds = pd.DataFrame(temp)
    if df_cmds.shape[0] != df_cmds["sample"].nunique():
        scp.logg.error("Samples are not unique!")

    cmd1 = step1(afterok)
    code1 = subprocess.getoutput(cmd1)
    cmd2 = step2(code1)
    code2 = subprocess.getoutput(cmd2)
    cmd3 = step3(code2)
    code3 = subprocess.getoutput(cmd3)
    cmd4 = step4(code3)
    code4 = subprocess.getoutput(cmd4)
    cmd5 = step5(code4)
    code5 = subprocess.getoutput(cmd5)
    scp.logg.info(code5)

    return None
