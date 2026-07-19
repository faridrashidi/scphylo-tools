"""Exercise external solver adapters without launching their executables."""

from tempfile import TemporaryDirectory
from unittest.mock import Mock

import networkx as nx
import pandas as pd
import pytest

from scphylo.tl.solver import _cellphy, _gpps, _grmt, _sciphi, _scite, _siclonefit


def test_cellphy_writes_vcf_and_reads_mocked_tree(monkeypatch):
    """Cover CellPhy command orchestration and its small-input warning."""
    tmpdir = TemporaryDirectory()
    tree_path = f"{tmpdir.name}/cellphy.vcf.raxml.bestTree"
    with open(tree_path, "w") as stream:
        stream.write("(cell1,cell2);\n")
    adata = Mock()
    adata.shape = (2, 2)
    monkeypatch.setattr(_cellphy.scp.ul, "executable", Mock(return_value="cellphy.sh"))
    monkeypatch.setattr(_cellphy.scp.ul, "tmpdirsys", Mock(return_value=tmpdir))
    monkeypatch.setattr(_cellphy.scp.io, "to_vcf", Mock())
    monkeypatch.setattr(_cellphy.scp.logg, "error", Mock())
    monkeypatch.setattr(_cellphy.scp.logg, "info", Mock())
    monkeypatch.setattr(_cellphy.os, "system", Mock(return_value=0))
    monkeypatch.setattr(_cellphy.time, "time", Mock(side_effect=[1.0, 2.0]))

    assert _cellphy.cellphy(adata, mode="full") == "(cell1,cell2);"
    _cellphy.scp.logg.error.assert_called_once()
    assert "FULL" in _cellphy.os.system.call_args.args[0]


def test_grmt_converts_mocked_dot_tree_to_matrix(monkeypatch):
    """Cover GRMT input staging and mutation-path conversion."""
    tmpdir = TemporaryDirectory()
    matrix = pd.DataFrame([[3], [3]], index=["cell1", "cell2"], columns=["mut1"])
    tree = nx.DiGraph()
    tree.add_node("0", label='"root"')
    tree.add_node("1", label='"mut1"')
    tree.add_node("2", label='"cell1"')
    tree.add_node("3", label='"cell2"')
    tree.add_edges_from([("0", "1"), ("1", "2"), ("0", "3")])
    monkeypatch.setattr(_grmt.scp.ul, "executable", Mock(return_value="grmt"))
    monkeypatch.setattr(_grmt.scp.ul, "tmpdirsys", Mock(return_value=tmpdir))
    monkeypatch.setattr(_grmt.scp.ul, "stat", Mock())
    monkeypatch.setattr(_grmt.scp.logg, "info", Mock())
    monkeypatch.setattr(_grmt.os, "system", Mock(return_value=0))
    monkeypatch.setattr(_grmt.time, "time", Mock(side_effect=[1.0, 2.0]))
    monkeypatch.setattr(_grmt.nx.drawing.nx_pydot, "read_dot", Mock(return_value=tree))

    result = _grmt.grmt(matrix, 0.1, 0.2, n_iters=2, n_threads=3)

    assert result.loc["cell1", "mut1"] == 1
    assert result.loc["cell2", "mut1"] == 0
    assert "--n_iter 2" in _grmt.os.system.call_args.args[0]


def test_sciphi_stages_all_genotype_symbols(monkeypatch, tmp_path):
    """Build SCIPhI's mpileup representation in an isolated directory."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "test").mkdir()
    matrix = pd.DataFrame(
        [[0, 1, 3], [1, 3, 0]],
        index=["cell1", "cell2"],
        columns=["m1", "m2", "m3"],
    )
    monkeypatch.setattr(_sciphi.scp.ul, "executable", Mock(return_value="sciphi"))
    monkeypatch.setattr(_sciphi.scp.logg, "info", Mock())
    monkeypatch.setattr(_sciphi.os, "system", Mock(return_value=0))
    monkeypatch.setattr(_sciphi.time, "time", Mock(side_effect=[1.0, 2.0]))

    assert _sciphi.sciphi(matrix) is None
    pileup = (tmp_path / "test" / "sciphi.mpileup").read_text()
    assert "\t.\t<" in pileup
    assert "\tT\t<" in pileup
    assert "\tN\t<" in pileup
    assert (tmp_path / "test" / "sciphi.cellnames").read_text().splitlines() == [
        "cell1\tCT",
        "cell2\tCT",
    ]


def test_gpps_swaps_error_rates_and_preserves_labels(monkeypatch):
    """Cover GPPS composition while replacing ILP and hill climbing."""
    matrix = pd.DataFrame(
        [[0, 1], [1, 0]], index=["cell1", "cell2"], columns=["m1", "m2"]
    )
    ilp = Mock(return_value=[[0, 1], [1, 0]])
    hill_climb = Mock(return_value=[[0, 1], [1, 0]])
    monkeypatch.setattr(_gpps, "gpps_ilp", ilp)
    monkeypatch.setattr(_gpps, "gpps_hc", hill_climb)
    monkeypatch.setattr(_gpps.time, "time", Mock(side_effect=[1.0, 2.5]))
    monkeypatch.setattr(_gpps.scp.ul, "stat", Mock())
    monkeypatch.setattr(_gpps.scp.logg, "info", Mock())

    result = _gpps.gpps(
        matrix,
        0.1,
        0.2,
        k_dollo=1,
        max_del=2,
        neighbor_size=3,
        n_iters=4,
        time_limit=5,
        n_threads=6,
        tree_script="tree.rb",
        ruby_executable="ruby",
        tools_dir="tools",
    )

    assert result.equals(matrix.rename_axis("cellIDxmutID"))
    assert ilp.call_args.kwargs["alpha"] == 0.2
    assert ilp.call_args.kwargs["beta"] == 0.1
    assert hill_climb.call_args.kwargs["tree_script"] == "tree.rb"


def test_infscite_converts_mocked_outputs_in_both_modes(monkeypatch):
    """Replace the former placeholder with a complete in-process infSCITE test."""
    matrix = pd.DataFrame([[0], [1]], index=["cell1", "cell2"], columns=["m1"])
    tempdirs = [TemporaryDirectory(), TemporaryDirectory()]
    for tmpdir in tempdirs:
        with open(f"{tmpdir.name}/infscite_ml0.gv", "w") as stream:
            stream.write("2 -> m1;\nm1 -> 1;\n")
        with open(f"{tmpdir.name}/infscite.log", "w") as stream:
            stream.write(
                "best value for beta: 0.2\n"
                "best value for alpha: 0.1\n"
                "best doublet rate: 0.05\n"
                "best log score for tree: -3.0\n"
            )
    tree = nx.DiGraph([("2", "m1"), ("m1", "1")])
    monkeypatch.setattr(_scite.scp.ul, "tmpdirsys", Mock(side_effect=tempdirs))
    monkeypatch.setattr(_scite.scp.ul, "get_file", Mock(return_value="infSCITE"))
    monkeypatch.setattr(_scite.scp.ul, "stat", Mock())
    monkeypatch.setattr(_scite.scp.logg, "info", Mock())
    monkeypatch.setattr(_scite.os, "system", Mock(return_value=0))
    monkeypatch.setattr(_scite.time, "time", Mock(side_effect=[1.0, 2.0, 3.0, 5.0]))
    monkeypatch.setattr(_scite.nx.drawing.nx_pydot, "read_dot", Mock(return_value=tree))

    regular = _scite.infscite(matrix, 0.1, 0.2, 2, n_restarts=1)
    experiment = _scite.infscite(matrix, 0.1, 0.2, 2, n_restarts=1, experiment=True)

    assert regular.loc["cell1", "m1"] == 1
    assert regular.loc["cell2", "m1"] == 0
    assert experiment[1:] == (2.0, -3.0, 0.2)


def test_siclonefit_translates_output_failures_and_experiment(monkeypatch, tmp_path):
    """Cover missing, unreadable, malformed, tree, and experiment result paths."""
    matrix = pd.DataFrame([[0, 1], [1, 0]], index=["c1", "c2"], columns=["a", "b"])
    jar = tmp_path / "SiCloneFiTComplete.jar"
    java = tmp_path / "java"
    jar.write_text("fixture")
    java.write_text("#!/bin/sh\n")
    java.chmod(0o755)
    modes = iter(["missing", "unreadable", "wrong-shape", "missing-tree", "valid"])

    def fake_run(args, _appname, **_kwargs):
        mode = next(modes)
        if mode == "missing":
            return None
        workdir = args[args.index("-outDir") + 1]
        best = workdir / "fixture_samples" / "best"
        best.mkdir(parents=True)
        genotype = best / "best_MAP_predicted_genotype.txt"
        if mode == "unreadable":
            genotype.mkdir()
        elif mode == "wrong-shape":
            pd.DataFrame([[0]]).to_csv(genotype, sep=" ", header=False)
        else:
            output = matrix.T.copy()
            output.index = range(output.shape[0])
            output.to_csv(genotype, sep=" ", header=False)
        return None

    monkeypatch.setattr(_siclonefit.scp.ul, "run_external", fake_run)
    monkeypatch.setattr(_siclonefit.time, "perf_counter", Mock(side_effect=range(10)))
    monkeypatch.setattr(_siclonefit.scp.ul, "stat", Mock())
    monkeypatch.setattr(
        _siclonefit.scp.ul, "is_conflict_free_gusfield", Mock(return_value=True)
    )
    monkeypatch.setattr(_siclonefit.scp.ul, "calc_nll_matrix", Mock(return_value=1.5))

    common = {
        "alpha": 0.1,
        "beta": 0.2,
        "jar_path": jar,
        "java_executable": java,
    }
    for message in (
        "without producing",
        "unreadable genotype",
        "shape",
        "readable MAP tree",
    ):
        with pytest.raises(
            _siclonefit.scp.ul.ExternalToolExecutionError, match=message
        ):
            _siclonefit.siclonefit(matrix, return_tree=True, **common)

    output, running_time, is_cf, nll = _siclonefit.siclonefit(
        matrix, experiment=True, **common
    )
    pd.testing.assert_frame_equal(output, matrix.rename_axis("cellIDxmutID"))
    assert (running_time, is_cf, nll) == (1, True, 1.5)
