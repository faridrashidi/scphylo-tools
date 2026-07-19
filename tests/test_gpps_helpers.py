"""Fast tests for GPPS tree helpers and its solver adapter."""

from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pandas as pd
import pytest

from scphylo.external.gpps import _gpps_hc, _gpps_ilp, _nh2lgf, _utils_hc


def _small_tree():
    """Return a tiny gain/loss tree and node lookup."""
    root = _utils_hc.Node("germline", None, 0, -1, tot_mutations=2)
    gain = _utils_hc.Node("m1", root, 1, 0)
    second = _utils_hc.Node("m2", gain, 2, 1)
    loss = _utils_hc.Node("m1---", gain, 3, 0, loss=True)
    return root, {node.id_node: node for node in (root, gain, second, loss)}


def test_newick_parser_handles_nested_lengths_and_invalid_input():
    """Exercise leaf and recursive parsing, with and without branch lengths."""
    nodes, edges = _nh2lgf.newick_to_edgelist("((m2:2)m1:1)m0;")
    assert set(nodes.values()) == {"m0", "m1", "m2"}
    assert len(edges) == 2
    assert _nh2lgf.newick_to_edgelist("leaf;") == ({0: "leaf"}, [])

    with pytest.raises(AssertionError):
        _nh2lgf.newick_to_edgelist("(broken;")
    with pytest.raises(AssertionError):
        _nh2lgf.node("(leaf);", 0, [], [])


def test_tree_nodes_copy_relationships_and_loss_cleanup(tmp_path, capsys):
    """Cover genotype construction, copying, movement, and invalid loss removal."""
    root, nodes = _small_tree()
    gain, second, loss = nodes[1], nodes[2], nodes[3]
    assert gain.is_ancestor_of(second)
    assert not second.is_ancestor_of(gain)
    assert _utils_hc.contains(gain.genotype_profile, root.genotype_profile)
    assert not _utils_hc.contains(root.genotype_profile, gain.genotype_profile)
    assert _utils_hc.is_loss_valid(loss, 0)
    assert not _utils_hc.is_already_lost(loss, 0)

    copied, copied_nodes = _utils_hc.copy_tree(root)
    assert copied is not root
    assert copied_nodes[2].genotype_profile == second.genotype_profile
    assert not _utils_hc.prune_and_reattach(
        copied_nodes[1], copied_nodes[2], copied_nodes
    )
    assert _utils_hc.prune_and_reattach(copied_nodes[2], copied_nodes[0], copied_nodes)

    invalid_loss = _utils_hc.Node("m2---", root, 4, 1, loss=True)
    nodes[4] = invalid_loss
    _utils_hc.check_subtree_losses(invalid_loss, nodes)
    assert 4 not in nodes

    repeated = _utils_hc.Node("m1---", loss, 5, 0, loss=True)
    nodes[5] = repeated
    assert _utils_hc.is_already_lost(repeated, 0)
    _utils_hc.check_subtree_losses(repeated, nodes)
    assert 5 not in nodes

    output = tmp_path / "tree.dot"
    with output.open("w") as stream:
        _utils_hc.print_dot_tree_file(root, stream)
    _utils_hc.print_dot_tree(root)
    dot_file = output.read_text()
    dot_stdout = capsys.readouterr().out
    assert '"0" -> "1"' in dot_file
    assert "color=indianred1" in dot_file
    assert '"0" -> "1"' in dot_stdout
    assert "color=indianred1" in dot_stdout


def test_build_and_import_tree_use_resolved_external_tool(monkeypatch):
    """Build a mutation tree from a mocked Ruby Newick conversion."""
    matrix = pd.DataFrame([[1, 0]])
    monkeypatch.setattr(
        _utils_hc.scp.ul, "resolve_external_file", Mock(return_value="tree.rb")
    )
    monkeypatch.setattr(
        _utils_hc.scp.ul, "resolve_executable", Mock(return_value="ruby")
    )
    run = Mock(return_value=SimpleNamespace(stdout="(m1,m2)germline;"))
    monkeypatch.setattr(_utils_hc.scp.ul, "run_external", run)

    root, nodes = _utils_hc.build_tree_from_file(matrix, ["gain", "loss---"], [0, 0], 1)
    assert root.name == "germline"
    assert len(nodes) == 3
    assert nodes[0].genotype_profile in ([1], [0])
    run.assert_called_once()

    monkeypatch.setattr(
        _utils_hc.scp.ul,
        "run_external",
        Mock(return_value=SimpleNamespace(stdout="")),
    )
    with pytest.raises(_utils_hc.scp.ul.ExternalToolExecutionError):
        _utils_hc.build_tree_from_file(matrix, ["gain", "loss---"], [0, 0], 1)

    build = Mock(return_value=(root, nodes))
    monkeypatch.setattr(_utils_hc, "build_tree_from_file", build)
    assert _utils_hc.import_ilp_out(matrix, 1, ["A"])[0] is root
    assert build.call_args.args[1:4] == (["A", "A---"], [0, 0], 1)


def test_hill_climbing_helpers_score_attach_and_improve(monkeypatch):
    """Exercise likelihood branches, attachments, neighbors, and improvement."""
    _gpps_hc.cell_row_likelihood.cache_clear()
    assert np.isfinite(_gpps_hc.cell_row_likelihood("012", "012", 0.1, 0.2))
    assert np.isfinite(_gpps_hc.cell_row_likelihood("0", "1", 0.1, 0.2))
    assert np.isfinite(_gpps_hc.cell_row_likelihood("1", "0", 0.1, 0.2))
    assert _gpps_hc.cell_row_likelihood("0", "2", 0.1, 0.2) == -np.inf
    assert _gpps_hc.cell_row_likelihood("1", "2", 0.1, 0.2) == -np.inf

    root, nodes = _small_tree()
    matrix = np.array([[0, 0], [1, 0], [1, 1], [2, 2]])
    likelihood, attachments = _gpps_hc.greedy_tree_likelihood(
        root, nodes, matrix, 0.1, 0.2
    )
    assert np.isfinite(likelihood)
    assert len(attachments) == len(matrix)
    expected = _gpps_hc.get_expect_matrix(root, nodes, matrix, 0.1, 0.2)
    assert len(expected) == len(matrix)

    choices = iter([2, 0])
    monkeypatch.setattr(_gpps_hc.random, "choice", lambda _items: next(choices))
    neighbors = _gpps_hc.generate_neighborhood(root, nodes, 1)
    assert len(neighbors) == 1

    better_root = _utils_hc.Node("germline", None, 10, -1, tot_mutations=2)
    better = _utils_hc.Node("perfect", better_root, 11, 0)
    better.genotype_profile = [1, 1]
    better_nodes = {10: better_root, 11: better}
    monkeypatch.setattr(
        _gpps_hc,
        "generate_neighborhood",
        Mock(return_value=[(better_root, better_nodes)]),
    )
    result = _gpps_hc.hill_climbing(
        root, {root.id_node: root}, 1, 2, 0.1, 0.2, np.array([[1, 1]])
    )
    assert result == (better_root, better_nodes)


def test_gpps_hc_composes_import_search_and_expected_matrix(monkeypatch):
    """Cover the high-level hill-climbing composition with tiny fakes."""
    root, nodes = _small_tree()
    monkeypatch.setattr(_gpps_hc, "import_ilp_out", Mock(return_value=(root, nodes)))
    monkeypatch.setattr(_gpps_hc, "hill_climbing", Mock(return_value=(root, nodes)))
    expected = [[0, 0]]
    monkeypatch.setattr(_gpps_hc, "get_expect_matrix", Mock(return_value=expected))

    assert (
        _gpps_hc.gpps_hc(
            np.array([[0, 0]]),
            pd.DataFrame([[0, 0]]),
            0.1,
            0.2,
            1,
            ["A", "B"],
            ns=1,
            mi=1,
        )
        == expected
    )


class _Expression:
    """Minimal symbolic expression used by the fake Gurobi model."""

    def __init__(self, value=1):
        self.X = value

    def _result(self, *_args):
        return _Expression()

    __add__ = _result
    __radd__ = _result
    __sub__ = _result
    __rsub__ = _result
    __eq__ = _result
    __ge__ = _result
    __le__ = _result


class _Model:
    """Small in-memory stand-in for the Gurobi API used by GPPS."""

    def __init__(self, _name):
        self.Params = SimpleNamespace()
        self.modelSense = None
        self.constraints = []

    def setParam(self, *_args):
        return None

    def addVar(self, **_kwargs):
        return _Expression()

    def addConstr(self, expression, *_args, **_kwargs):
        self.constraints.append(expression)
        return expression

    def update(self):
        return None

    def optimize(self):
        return None


def test_gpps_ilp_builds_model_with_fake_gurobi(monkeypatch):
    """Exercise every variable/constraint family without a Gurobi license."""
    fake_gp = SimpleNamespace(
        Model=_Model,
        GRB=SimpleNamespace(BINARY=1, MAXIMIZE=2),
        quicksum=lambda values: sum(values, _Expression(0)),
    )
    monkeypatch.setattr(
        _gpps_ilp.scp.ul, "import_gurobi", Mock(return_value=(fake_gp, True))
    )
    monkeypatch.setattr(_gpps_ilp.scp.logg, "error", Mock())
    monkeypatch.setattr(_gpps_ilp.scp.logg, "debug", Mock())

    result = _gpps_ilp.gpps_ilp(
        [[0, 1, 2], [1, 0, 1]],
        alpha=0.1,
        beta=0.2,
        k_dollo=1,
        max_del=1,
        time_limit=3,
        n_threads=2,
    )

    assert result == [[1] * 6, [1] * 6]
    _gpps_ilp.scp.logg.error.assert_called_once()
