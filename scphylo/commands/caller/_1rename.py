import glob
import os

import click

import scphylo as scp


@click.command(short_help="Run Rename.")
@click.argument(
    "indir",
    required=True,
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, readable=True, resolve_path=True
    ),
)
@click.argument(
    "keep",
    required=True,
    type=str,
)
@click.option(
    "--exec",
    is_flag=True,
    default=False,
    type=bool,
    show_default=True,
    help="Execute the renaming.",
)
def rename(indir, keep, exec):
    """Run Rename.

    trisicell rename path/to/in/dir 0,1,3 --exec
    """
    keep = list(map(int, keep.split(",")))
    idx = 0
    for file1 in sorted(glob.glob(indir + "/*.fastq.gz")):
        dirname, basename = scp.ul.dir_base(file1)
        file2 = [x for i, x in enumerate(basename.split("_")) if i in keep]
        file2 = "_".join(file2)
        if idx % 2 == 0:
            file2 += "_1.fastq.gz"
        else:
            file2 += "_2.fastq.gz"
        file2 = dirname + "/" + file2
        if not exec:
            scp.logg.info(file1, file2)
        else:
            os.rename(file1, file2)
        idx += 1
    return None
