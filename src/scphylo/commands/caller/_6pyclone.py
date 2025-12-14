import glob
import subprocess

import click
import pandas as pd

import scphylo as scp
from scphylo.ul._servers import cmd, write_cmds_get_main


@click.command(short_help="Run PyClone.")
@click.argument(
    "outdir",
    required=True,
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, readable=True, resolve_path=True
    ),
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
def pyclone(outdir, time, mem, afterok):
    """Run PyClone.

    scphylo caller pyclone /path/to/sequenza/folders

    Sequenza folders (*.sequenza) --> PyClone files
    """

    def get_command(sample):
        cmds = ""
        cmds += cmd(["unset PATH"])
        cmds += cmd(["unset PYTHONPATH"])
        cmds += cmd([f"module load {scp.settings.tools['pyclone']}"])
        cmds += cmd(["module load python/2.7"])
        with open(f"{outdir}/{sample}.sequenza/{sample}_cellularity.txt") as file:
            tumor_content = file.readline().strip()
        cmds += cmd(
            [
                "PyClone setup_analysis",
                f"--in_files {outdir}/{sample}.sequenza/{sample}_pyclone.tsv",
                f"--tumour_contents {tumor_content}",
                "--prior major_copy_number",
                f"--working_dir {outdir}/{sample}.pyclone",
            ]
        )
        cmds += cmd(
            [
                "PyClone run_analysis",
                f"--config_file {outdir}/{sample}.pyclone/config.yaml",
            ]
        )
        cmds += cmd(
            [
                "PyClone build_table",
                f"--config_file {outdir}/{sample}.pyclone/config.yaml",
                f"--out_file {outdir}/{sample}.pyclone/loci.tsv",
                "--table_type loci",
            ]
        )
        cmds += cmd(
            [
                "PyClone plot_clusters",
                f"--config_file {outdir}/{sample}.pyclone/config.yaml",
                f"--plot_file {outdir}/{sample}.pyclone/coordinates_cluster_plot.pdf",
                "--plot_type parallel_coordinates",
            ]
        )
        cmds += cmd(
            [
                "PyClone plot_clusters",
                f"--config_file {outdir}/{sample}.pyclone/config.yaml",
                f"--plot_file {outdir}/{sample}.pyclone/density_cluster_plot.pdf",
                "--plot_type density",
            ]
        )
        cmds += cmd(["echo Done!"], islast=True)
        return cmds

    files = glob.glob(f"{outdir}/*.sequenza")
    temp = []
    for file in files:
        file = file.split("/")[-1].replace(".sequenza", "")
        temp.append({"sample": file})
    df_cmds = pd.DataFrame(temp)
    df_cmds["cmd"] = df_cmds.apply(lambda x: get_command(x["sample"]), axis=1)

    cmdmain = write_cmds_get_main(
        df_cmds,
        "pyclone",
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
