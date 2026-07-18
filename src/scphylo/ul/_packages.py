import scphylo as scp


def import_gurobi():
    try:
        import gurobipy as gp

        return gp, False
    except ImportError:
        scp.logg.warn(
            "Unable to import `gurobipy`!",
            "Make sure `Gurobi` is already installed in your system.",
            "Then install `gurobipy`.",
        )
        return None, True


def import_mpi4py():
    try:
        import mpi4py

        return mpi4py, False
    except ImportError:
        scp.logg.warn(
            "Unable to import `mpi4py`!",
            "Install the locked Pixi development environment, or install an MPI "
            "implementation and `mpi4py` manually.",
        )
        return None, True


def import_rpy2(name="base", how=""):
    try:
        import logging

        from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
        from rpy2.robjects import r
        from rpy2.robjects.packages import PackageNotInstalledError, importr

        type(r)
        rpy2_logger.setLevel(logging.ERROR)

    except ImportError:
        scp.logg.warn(
            "Unable to import `rpy2`! Install the locked Pixi plotting environment "
            "with `pixi install`, or install R and `rpy2` manually."
        )
        return None, True

    try:
        _r_lib = importr(name)
        return _r_lib, False
    except PackageNotInstalledError:
        scp.logg.warn(f"Install R library `{name!r}` first.\n{how}")
        return None, True


def import_graphviz():
    try:
        import pygraphviz

        return pygraphviz, False
    except ImportError:
        scp.logg.warn(
            "Unable to import `pygraphviz`!",
            "Install the locked Pixi plotting environment with `pixi install`, or "
            "install both native Graphviz and `pygraphviz` manually.",
        )
        return None, True


def import_graph_tool():
    try:
        import graph_tool

        return graph_tool, False
    except ImportError:
        scp.logg.warn(
            "Unable to import `graph_tool`!",
            "Install the locked Pixi development environment.",
        )
        return None, True
