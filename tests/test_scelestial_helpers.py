"""Cover Scelestial's pure-Python format conversion and orchestration."""

from tempfile import TemporaryDirectory
from unittest.mock import Mock

import graphviz
import networkx as nx
import pandas as pd
import pytest
from scipy import stats

from scphylo.tl.solver import _scelestial as scelestial_module


def test_scelestial_orchestrates_converters_without_binary(monkeypatch):
    """Run the public adapter with all external work replaced by tiny fakes."""
    tmpdir = TemporaryDirectory()
    input_matrix = pd.DataFrame(
        [[0, 1], [1, 0]], index=["cell1", "cell2"], columns=["mut1", "mut2"]
    )
    dot_tree = nx.DiGraph()
    dot_tree.add_node("0", label='""')
    dot_tree.add_node("1", label="cell")
    dot_tree.add_edge("0", "1")
    output_matrix = pd.DataFrame([[0], [1]], index=input_matrix.index, columns=["m"])
    monkeypatch.setattr(
        scelestial_module.scp.ul, "executable", Mock(return_value="scelestial")
    )
    monkeypatch.setattr(
        scelestial_module.scp.ul, "tmpdirsys", Mock(return_value=tmpdir)
    )
    monkeypatch.setattr(scelestial_module.scp.ul, "root_id", Mock(return_value=0))
    monkeypatch.setattr(
        scelestial_module.scp.ul, "to_cfmatrix", Mock(return_value=output_matrix)
    )
    monkeypatch.setattr(scelestial_module.scp.ul, "stat", Mock())
    monkeypatch.setattr(scelestial_module.scp.logg, "info", Mock())
    monkeypatch.setattr(scelestial_module.os, "system", Mock(return_value=0))
    monkeypatch.setattr(scelestial_module.time, "time", Mock(side_effect=[1.0, 2.5]))
    monkeypatch.setattr(scelestial_module, "_convert_input", Mock())
    monkeypatch.setattr(scelestial_module, "_steiner_to_seq", Mock())
    monkeypatch.setattr(scelestial_module, "_stein_to_clone_tree", Mock())
    monkeypatch.setattr(scelestial_module, "_clone_tree_to_mu_tree_imput", Mock())
    monkeypatch.setattr(
        scelestial_module.nx.nx_pydot, "read_dot", Mock(return_value=dot_tree)
    )

    result = scelestial_module.scelestial(input_matrix)

    assert result.shape == (2, 3)
    assert result.columns.tolist() == ["m", "mut2", "mut3"]
    scelestial_module.scp.ul.stat.assert_called_once_with(result, result, 0, 0, 1.5)


def test_convert_input_handles_every_genotype_and_rejects_unknown(tmp_path, capsys):
    """Convert all supported genotype symbols and cover validation."""
    source = tmp_path / "input.SC.T"
    imputed = tmp_path / "input.txt"
    breakpoints = tmp_path / "breakpoints.txt"
    source.write_text("0 1 2 3\n1 0 3 2\n")

    scelestial_module._convert_input(source, imputed, breakpoints)

    assert "A/A" in imputed.read_text()
    assert "C/C" in imputed.read_text()
    assert "A/C" in imputed.read_text()
    assert "./." in imputed.read_text()
    assert breakpoints.read_text().splitlines()[0] == "V1,V2"

    source.write_text("4\n")
    with pytest.raises(ValueError, match="OH! 4"):
        scelestial_module._convert_input(source, imputed, breakpoints)
    assert "OH! 4" in capsys.readouterr().out


def test_steiner_converters_cover_tree_and_imputation_paths(
    monkeypatch, tmp_path, capsys
):
    """Exercise clone-tree rooting and valid, empty, and invalid imputations."""
    steiner = tmp_path / "steiner.txt"
    steiner.write_text("3\n0 1 AAA ACA\n1 0 ACC ACC\n2 1 CCC CCA\n2\n0 1 1\n1 2 2\n")
    sequence = tmp_path / "sequence.txt"
    sequence.write_text("unused\n")
    tree = tmp_path / "tree.txt"
    clones = tmp_path / "clones.txt"
    monkeypatch.setattr(
        scelestial_module.sys,
        "argv",
        ["scelestial", "1", "2", "3", "4", "-exclude-root"],
    )

    scelestial_module._stein_to_clone_tree(sequence, steiner, tree, clones)
    assert tree.read_text().splitlines()[0] == "0 1 2"
    assert clones.read_text().splitlines()[0] == "0 "

    imputed = tmp_path / "imputed.txt"
    scelestial_module._steiner_to_seq(steiner, imputed)
    assert imputed.read_text().splitlines() == ["0 1", "1 1", "0 0"]

    steiner.write_text("1\n0 0 AAA AAA\n")
    with pytest.raises(SystemExit):
        scelestial_module._steiner_to_seq(steiner, imputed)

    steiner.write_text("1\n0 1 AAA AXA\n")
    with pytest.raises(Exception, match="Invalid imputation"):
        scelestial_module._steiner_to_seq(steiner, imputed)
    assert "invalid char" in capsys.readouterr().err


def test_clone_tree_rendering_options_use_graphviz_fake(monkeypatch, tmp_path):
    """Cover mutation labels, compression, separated plots, and clone merging."""
    tree = tmp_path / "tree.txt"
    clones = tmp_path / "clones.txt"
    sequences = tmp_path / "sequences.txt"
    mutations = tmp_path / "mutations.txt"
    cells = tmp_path / "cells.txt"
    tree.write_text("0 1 2\n1->0 1\n2->1 1\n")
    clones.write_text("0 1\n1 2\n2 3\n")
    sequences.write_text("0 1 1\n0 0 1\n")
    mutations.write_text("mut1\nmut2\n")
    cells.write_text("cell1\ncell2\ncell3\n")

    renders = []

    class FakeDigraph:
        """Record nodes, edges, and renders without invoking Graphviz."""

        def __init__(self, **_kwargs):
            self.graph_attr = {}

        def node(self, *_args, **_kwargs):
            return None

        def edge(self, *_args, **_kwargs):
            return None

        def render(self, output):
            renders.append(output)

    monkeypatch.setattr(graphviz, "Digraph", FakeDigraph)
    monkeypatch.setattr(stats, "binom_test", Mock(return_value=0.0), raising=False)

    scelestial_module._clone_tree_to_mu_tree_imput(
        tree,
        clones,
        sequences,
        mutations,
        cells,
        str(tmp_path / "output"),
        markMutations=True,
        compress=True,
        markMutationsSeparated=True,
    )
    scelestial_module._clone_tree_to_mu_tree_imput(
        tree,
        clones,
        sequences,
        mutations,
        cells,
        str(tmp_path / "merged"),
        margeClones=True,
    )

    assert renders == [
        str(tmp_path / "output-0"),
        str(tmp_path / "output-1"),
        str(tmp_path / "merged"),
    ]
