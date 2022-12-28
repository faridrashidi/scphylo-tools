import pytest

import scphylo as scp

skip_gurobi = pytest.mark.skipif(
    scp.ul.import_gurobi()[1], reason="Unable to import `Gurobi`!"
)
skip_mpi4py = pytest.mark.skipif(
    scp.ul.import_mpi4py()[1], reason="Unable to import `MPI`!"
)
skip_graphviz = pytest.mark.skipif(
    scp.ul.import_graphviz()[1], reason="Unable to import `Graphviz`!"
)
skip_graph_tool = pytest.mark.skipif(
    scp.ul.import_graph_tool()[1], reason="Unable to import `graph_tool`!"
)


def skip_rpy2(package="base"):
    return pytest.mark.skipif(
        scp.ul.import_rpy2(package)[1], reason="Unable to import `rpy2`!"
    )
