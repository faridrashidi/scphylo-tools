"""Cover lightweight MP3 and hierarchical-clustering edge cases."""

from collections import Counter

import networkx as nx
import numpy as np
import pandas as pd
import pytest

from scphylo.external import _mp3
from scphylo.ul import _hclustering


def _star_tree(labels):
    """Return an MP3 tree whose requested labels are sibling leaves."""
    graph = nx.DiGraph()
    graph.add_node(("root",))
    for label in labels:
        node = (label,)
        graph.add_edge(("root",), node)
        graph.nodes[node]["label"] = label
    return _mp3.build_tree(graph, labeled_only=True)


def test_mp3_lca_and_scalar_helpers():
    """Exercise direct LCA lookup and the bounded score helpers."""
    tree = _star_tree("ab")

    assert tree.LCA.lca_labels("a", "b") == Counter({"root": 1})
    assert "('root',)" in str(tree.LCA)
    assert _mp3.sigmoid(0) == 0
    assert _mp3.sigmoid(1) == 1
    assert _mp3.get_nset_sig(0, 0.25) == 0.25


def test_mp3_similarity_modes_and_validation():
    """Cover shared-label, union, geometric, and invalid score modes."""
    three_labels = _star_tree("abc")
    four_labels = _star_tree("abcd")

    with pytest.raises(AttributeError, match="Incorrect value"):
        _mp3.similarity(three_labels, four_labels, mode="invalid")

    assert _mp3.similarity(_star_tree("a"), _star_tree("z")) == 0.0
    assert _mp3.similarity(three_labels, four_labels, mode="intersection") == 1.0
    assert _mp3.similarity(three_labels, four_labels, mode="union") == 0.25
    assert _mp3.similarity(three_labels, four_labels, mode="geometric") == 0.5


def test_mp3_parallel_path_without_processes(monkeypatch, capsys):
    """Exercise result aggregation through a synchronous pool stand-in."""

    class InlinePool:
        """Implement the small context-manager surface used by similarity."""

        def __init__(self, processes):
            assert processes is None

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def map(self, function, values):
            return [function(value) for value in values]

    def already_configured(_method):
        raise RuntimeError

    monkeypatch.setattr(_mp3, "Pool", InlinePool)
    monkeypatch.setattr(_mp3, "set_start_method", already_configured)

    score = _mp3.similarity(
        _star_tree("abc"), _star_tree("abcd"), mode="union", cores=0
    )

    assert score == 0.25
    assert capsys.readouterr().out == "ERROR\n"


def test_mp3_build_tree_label_options():
    """Validate tree rejection, generated labels, omissions, and exclusions."""
    with pytest.raises(ValueError, match="Not a valid tree"):
        _mp3.build_tree(nx.DiGraph([(0, 1), (1, 0)]))

    generated = nx.DiGraph([(0, 1)])
    generated.nodes[1]["label"] = "keep,drop"
    tree = _mp3.build_tree(generated, exclude={"drop"})
    assert tree.label_set == ["0", "keep"]

    labeled_only = nx.DiGraph([("root", "leaf")])
    labeled_only.nodes["leaf"]["label"] = "leaf"
    tree = _mp3.build_tree(labeled_only, labeled_only=True)
    assert tree.label_set == ["leaf"]
    assert tree.T.nodes["root"]["label"] == ""


def test_distance_kernel_edge_cases_without_jit():
    """Run Numba kernels as Python to cover missing-data branches cheaply."""
    l1 = _hclustering._l1_ignore_na.py_func
    cosine = _hclustering._cosine_ignore_na.py_func

    assert l1(np.array([0.0, 3.0, np.nan]), np.array([1.0, 0.0, 1.0])) == 1
    assert np.isnan(l1(np.array([3.0, np.nan]), np.array([3.0, 1.0])))
    assert cosine(np.array([np.nan, 0.0]), np.array([1.0, 0.0])) == 1
    assert cosine(np.ones(3), np.ones(3)) == 0


def test_hclustering_return_distance_and_invalid_metric(monkeypatch):
    """Return a precomputed distance matrix and reject unknown metrics."""
    frame = pd.DataFrame([[0, 1], [1, 0]], index=["a", "b"])
    expected = np.array([[0.0, 2.0], [2.0, 0.0]])
    monkeypatch.setattr(_hclustering, "dist_l1_ignore_na", lambda _values: expected)

    assert _hclustering.hclustering(frame, return_dist=True) is expected
    with pytest.raises(RuntimeError):
        _hclustering.hclustering(frame, metric="unknown")
