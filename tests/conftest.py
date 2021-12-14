import shutil

import pytest

import scphylo as scp


def pytest_sessionstart(session):
    scp.ul.mkdir("./.pytest_cache")
    shutil.copy2(
        scp.ul.get_file("scphylo.datasets/test/test.tsv"),
        "./.pytest_cache",
    )


# removes overly verbose and useless logging errors for rpy2
# see: https://github.com/pytest-dev/pytest/issues/5502#issuecomment-647157873
def pytest_sessionfinish(session, exitstatus):
    import logging

    loggers = [logging.getLogger()] + list(logging.Logger.manager.loggerDict.values())
    for logger in loggers:
        handlers = getattr(logger, "handlers", [])
        for handler in handlers:
            logger.removeHandler(handler)

    scp.ul.cleanup("./.pytest_cache")


@pytest.fixture(scope="session")
def test_dir():
    return "./.pytest_cache"


@pytest.fixture(scope="session")
def test_data():
    return "./.pytest_cache/test.tsv"


@pytest.fixture(scope="session")
def test_cf_data_1():
    return scp.ul.get_file("scphylo.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")


@pytest.fixture(scope="session")
def test_cf_data_2():
    return scp.ul.get_file("scphylo.datasets/test/fp_1-fn_0.1-na_0.bnb.CFMatrix")


@pytest.fixture(scope="session")
def test_consensus_biorxiv_fig3b():
    return scp.ul.get_file("scphylo.datasets/test/consensus/biorxiv.fig3b.CFMatrix")


@pytest.fixture(scope="session")
def test_consensus_biorxiv_figs18a():
    return scp.ul.get_file("scphylo.datasets/test/consensus/biorxiv.figs18a.CFMatrix")


@pytest.fixture(scope="session")
def test_consensus_recomb_fig1a():
    return scp.ul.get_file("scphylo.datasets/test/consensus/recomb.fig1a.CFMatrix")


@pytest.fixture(scope="session")
def test_consensus_recomb_fig1b():
    return scp.ul.get_file("scphylo.datasets/test/consensus/recomb.fig1b.CFMatrix")


@pytest.fixture(scope="session")
def test_consensus_biorxiv_fig4b():
    return scp.ul.get_file("scphylo.datasets/test/consensus/biorxiv.fig4b.CFMatrix")


@pytest.fixture(scope="session")
def test_consensus_biorxiv_fig4c():
    return scp.ul.get_file("scphylo.datasets/test/consensus/biorxiv.fig4c.CFMatrix")


@pytest.fixture(scope="session")
def test_consensus_biorxiv_fig3c():
    return scp.ul.get_file("scphylo.datasets/test/consensus/biorxiv.fig3c.CFMatrix")


@pytest.fixture(scope="session")
def test_newick_1():
    return scp.ul.get_file("scphylo.datasets/test/input/T00.nwk")


@pytest.fixture(scope="session")
def test_newick_2():
    return scp.ul.get_file("scphylo.datasets/test/input/T06.nwk")
