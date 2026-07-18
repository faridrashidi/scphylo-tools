"""Provide shared fixtures and lifecycle hooks for the test suite."""

import shutil

import pytest

import scphylo as scp


def pytest_sessionstart(session):
    """Copy the shared genotype fixture into the session cache."""
    scp.ul.mkdir("./.pytest_cache")
    shutil.copy2(
        scp.ul.get_file("scphylo.datasets/test/test.tsv"),
        "./.pytest_cache",
    )


# removes overly verbose and useless logging errors for rpy2
# see: https://github.com/pytest-dev/pytest/issues/5502#issuecomment-647157873
def pytest_sessionfinish(session, exitstatus):
    """Remove rpy2 logging handlers and clean the session cache."""
    import logging

    loggers = [logging.getLogger()] + list(logging.Logger.manager.loggerDict.values())
    for logger in loggers:
        handlers = getattr(logger, "handlers", [])
        for handler in handlers:
            logger.removeHandler(handler)

    scp.ul.cleanup("./.pytest_cache")


@pytest.fixture(scope="session")
def test_dir():
    """Return the path to the shared test cache directory."""
    return "./.pytest_cache"


@pytest.fixture(scope="session")
def test_data():
    """Return the path to the cached genotype test matrix."""
    return "./.pytest_cache/test.tsv"


@pytest.fixture(scope="session")
def test_cf_data_1():
    """Return the path to the ground-truth conflict-free matrix."""
    return scp.ul.get_file("scphylo.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")


@pytest.fixture(scope="session")
def test_cf_data_2():
    """Return the path to the inferred conflict-free matrix."""
    return scp.ul.get_file("scphylo.datasets/test/fp_1-fn_0.1-na_0.bnb.CFMatrix")


@pytest.fixture(scope="session")
def test_consensus_biorxiv_fig3b():
    """Return the first BioRxiv figure 3 consensus fixture."""
    return scp.ul.get_file("scphylo.datasets/test/consensus/biorxiv.fig3b.CFMatrix")


@pytest.fixture(scope="session")
def test_consensus_biorxiv_figs18a():
    """Return the BioRxiv supplementary figure 18 consensus fixture."""
    return scp.ul.get_file("scphylo.datasets/test/consensus/biorxiv.figs18a.CFMatrix")


@pytest.fixture(scope="session")
def test_consensus_recomb_fig1a():
    """Return the first RECOMB figure 1 consensus fixture."""
    return scp.ul.get_file("scphylo.datasets/test/consensus/recomb.fig1a.CFMatrix")


@pytest.fixture(scope="session")
def test_consensus_recomb_fig1b():
    """Return the second RECOMB figure 1 consensus fixture."""
    return scp.ul.get_file("scphylo.datasets/test/consensus/recomb.fig1b.CFMatrix")


@pytest.fixture(scope="session")
def test_consensus_biorxiv_fig4b():
    """Return the first BioRxiv figure 4 consensus fixture."""
    return scp.ul.get_file("scphylo.datasets/test/consensus/biorxiv.fig4b.CFMatrix")


@pytest.fixture(scope="session")
def test_consensus_biorxiv_fig4c():
    """Return the second BioRxiv figure 4 consensus fixture."""
    return scp.ul.get_file("scphylo.datasets/test/consensus/biorxiv.fig4c.CFMatrix")


@pytest.fixture(scope="session")
def test_consensus_biorxiv_fig3c():
    """Return the second BioRxiv figure 3 consensus fixture."""
    return scp.ul.get_file("scphylo.datasets/test/consensus/biorxiv.fig3c.CFMatrix")


@pytest.fixture(scope="session")
def test_newick_1():
    """Return the path to the first Newick tree fixture."""
    return scp.ul.get_file("scphylo.datasets/test/input/T00.nwk")


@pytest.fixture(scope="session")
def test_newick_2():
    """Return the path to the second Newick tree fixture."""
    return scp.ul.get_file("scphylo.datasets/test/input/T06.nwk")
