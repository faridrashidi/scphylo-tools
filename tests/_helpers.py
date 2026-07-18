import pytest

import scphylo as scp

skip_gurobi = pytest.mark.skipif(
    scp.ul.import_gurobi()[1], reason="Unable to import `Gurobi`!"
)


def skip_rpy2(package="base"):
    return pytest.mark.skipif(
        scp.ul.import_rpy2(package)[1], reason="Unable to import `rpy2`!"
    )
