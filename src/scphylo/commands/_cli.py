"""Command Line Interface."""

import click

import scphylo as scp


@click.group()
@click.version_option()
def cli():
    """scphylo-tools command line interface."""
    pass


@cli.group()
def solver():
    """Solver commands."""
    pass


@cli.group()
def utils():
    """Run utility commands."""
    pass


@solver.command()
@click.argument("data")
@click.argument("alpha", type=float)
@click.argument("beta", type=float)
def scistree(data, alpha, beta):
    """Run ScisTree solver."""
    df = scp.io.read(data)
    scp.tl.scistree(df, alpha=alpha, beta=beta)


@solver.command()
@click.argument("data")
@click.argument("alpha", type=float)
@click.argument("beta", type=float)
def huntress(data, alpha, beta):
    """Run HUNTRESS solver."""
    df = scp.io.read(data)
    scp.tl.huntress(df, alpha=alpha, beta=beta)


@solver.command()
@click.argument("data")
@click.argument("alpha", type=float)
@click.argument("beta", type=float)
@click.option("-r", "--n_restarts", default=3, type=int)
@click.option("-l", "--n_iters", default=1000, type=int)
@click.option("-e", "--experiment", is_flag=True)
@click.option("-t", "--n_threads", default=1, type=int)
@click.option("-s", "--seed", default=1, type=int)
def scite(data, alpha, beta, n_restarts, n_iters, experiment, n_threads, seed):
    """Run SCITE solver."""
    df = scp.io.read(data)
    scp.tl.scite(df, alpha=alpha, beta=beta, n_restarts=n_restarts, n_iters=n_iters)


@solver.command()
@click.argument("data")
@click.argument("alpha", type=float)
@click.argument("beta", type=float)
def phiscsb(data, alpha, beta):
    """Run PhISCS-B solver."""
    df = scp.io.read(data)
    scp.tl.phiscsb(df, alpha=alpha, beta=beta)


@solver.command()
@click.argument("data")
@click.argument("alpha", type=float)
@click.argument("beta", type=float)
@click.option("--solver", default="scite")
@click.option("--n_samples", default=100, type=int)
@click.option("--sample_size", default=15, type=int)
@click.option("--n_jobs", default=1, type=int)
@click.option("--n_iterations", default=10000, type=int)
def booster(data, alpha, beta, solver, n_samples, sample_size, n_jobs, n_iterations):
    """Run Booster solver."""
    df = scp.io.read(data)
    scp.tl.booster(
        df,
        alpha=alpha,
        beta=beta,
        solver=solver,
        n_samples=n_samples,
        sample_size=sample_size,
        n_jobs=n_jobs,
        n_iterations=n_iterations,
    )


@solver.command()
@click.argument("data")
@click.option("-b", "--bounding", default="simulated")
def bnb(data, bounding):
    """Run BnB solver."""
    df = scp.io.read(data)
    scp.tl.bnb(df, bounding=bounding)


@utils.command()
@click.argument("file1")
@click.argument("file2")
@click.argument("output")
def consensus(file1, file2, output):
    """Consensus tree."""
    sc1 = scp.io.read(file1)
    sc2 = scp.io.read(file2)
    tree = scp.tl.consensus(sc1, sc2)
    df = scp.ul.to_cfmatrix(tree)
    scp.io.write(df, output)


@utils.command()
@click.argument("file")
def cf2newick(file):
    """Convert CFMatrix to Newick."""
    pass


@utils.command()
@click.argument("file")
@click.option("-p", "--processes", default=1, type=int)
def search(file, processes):
    """Search."""
    pass


@utils.command()
@click.argument("file1")
@click.argument("file2")
def score(file1, file2):
    """Score."""
    pass


@utils.command()
@click.argument("file")
def cf2tree(file):
    """Convert CFMatrix to Tree."""
    pass


@utils.command()
@click.argument("data")
@click.argument("alpha", type=float)
@click.argument("beta", type=float)
@click.option("--n_threads", default=1, type=int)
@click.option("--n_samples", default=100, type=int)
def partf(data, alpha, beta, n_threads, n_samples):
    """Partition function."""
    df = scp.io.read(data)
    scp.tl.partition_function(df, alpha=alpha, beta=beta, n_samples=n_samples)
