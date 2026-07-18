import os

import pytest

import scphylo as scp

skip_gurobi = pytest.mark.skipif(
    scp.ul.import_gurobi()[1], reason="Unable to import `Gurobi`!"
)
skip_slow = pytest.mark.skipif(
    os.environ.get("SCPHYLO_RUN_SLOW_TESTS") != "1",
    reason="Set SCPHYLO_RUN_SLOW_TESTS=1 to run slow integration tests.",
)


def skip_rpy2(package="base"):
    """Skip a test when the requested R package is unavailable through rpy2."""
    return pytest.mark.skipif(
        scp.ul.import_rpy2(package)[1],
        reason=f"R package `{package}` is unavailable through `rpy2`.",
    )
