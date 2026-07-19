"""Cover small tree-tool branches without invoking expensive solver backends."""

from decimal import Decimal
from unittest.mock import Mock

import networkx as nx
import numpy as np
import pandas as pd
import pytest

import scphylo as scp
from scphylo.tl.consensus import _consensus as consensus_module
from scphylo.tl.partition_function import _clt_sampler, _partition_function, _pf
from scphylo.tl.score import _others, _ours


def _matrix(values, *, cells=("c1", "c2", "c3"), muts=("a", "b")):
    """Build a small labeled genotype matrix."""
    return pd.DataFrame(values, index=cells, columns=muts)


@pytest.mark.parametrize(
    ("score", "ground", "solution"),
    [
        (
            _ours.gs,
            _matrix([[1], [0]], cells=("c1", "c2"), muts=("a",)),
            _matrix([[1], [0]], cells=("c1", "c2"), muts=("b",)),
        ),
        (
            _ours.gs,
            _matrix([[1]], cells=("c1",), muts=("a",)),
            _matrix([[1]], cells=("c2",), muts=("a",)),
        ),
        *[
            (
                score,
                _matrix([[1], [0]], cells=("c1", "c2"), muts=("a",)),
                _matrix([[1], [0]], cells=("c1", "c2"), muts=("b",)),
            )
            for score in (
                _ours.ad,
                _ours.dl,
                _ours.cc,
                _ours.mltd,
                _ours.tpted,
                _others.mp3,
                _others.caset,
                _others.disc,
            )
        ],
        (
            _others.rf,
            _matrix([[1]], cells=("c1",), muts=("a",)),
            _matrix([[1]], cells=("c2",), muts=("a",)),
        ),
    ],
)
def test_scores_reject_inputs_without_shared_labels(score, ground, solution):
    """Reject trees that have no mutations or cells available for comparison."""
    with pytest.raises(RuntimeError):
        score(ground, solution)


@pytest.mark.parametrize(
    ("ground", "solution"),
    [
        (
            _matrix([[1, 1], [1, 0], [0, 0]]),
            _matrix([[1, 0], [0, 1], [0, 0]]),
        ),
        (
            _matrix([[1, 1], [0, 1], [0, 0]]),
            _matrix([[1, 1], [1, 0], [0, 0]]),
        ),
    ],
)
def test_ad_detects_each_reversed_ancestor_case(ground, solution):
    """Count disjoint and direction-reversed ancestor relationships as errors."""
    assert _ours.ad(ground, solution) == 0


def test_ad_rejects_a_tree_without_ancestor_pairs():
    """Reject a comparison for which the ground tree has no ancestor pair."""
    matrix = _matrix([[1], [0], [0]], muts=("a",))

    with pytest.raises(RuntimeError):
        _ours.ad(matrix, matrix)


def test_unimplemented_scores_return_none():
    """Document the current placeholder behavior of Bourque and PCSS scores."""
    matrix = _matrix([[1, 0], [0, 1], [0, 0]])

    assert _others.bourque(matrix, matrix) is None
    assert _others.pcss(matrix, matrix) is None


def test_mltd_reports_a_core_failure(monkeypatch):
    """Turn a failed native MLTD call into the package's standard runtime error."""
    matrix = _matrix([[1, 1], [1, 0], [0, 0]])
    monkeypatch.setattr(_ours, "run_mltd", Mock(return_value=None))

    with pytest.raises(RuntimeError):
        _ours.mltd(matrix.copy(), matrix.copy())


@pytest.mark.parametrize("function", [scp.tl.consensus, scp.tl.consensus_day])
def test_consensus_rejects_inputs_without_shared_nonzero_cells(function):
    """Reject consensus inputs without a shared tumor cell."""
    first = _matrix([[1]], cells=("c1",), muts=("a",))
    second = _matrix([[1]], cells=("c2",), muts=("b",))

    with pytest.raises(RuntimeError):
        function(first, second)


def test_consensus_contracts_two_labeled_nodes(monkeypatch):
    """Merge both cell labels when an unsupported internal edge is contracted."""
    first = nx.DiGraph(splitter_mut="|", splitter_cell="|")
    first.add_node(0, label="––")
    first.add_node(1, label="a")
    first.add_node(2, label="b")
    first.add_edge(0, 1, label="m0")
    first.add_edge(1, 2, label="m1")

    second = nx.DiGraph(splitter_mut="|", splitter_cell="|")
    second.add_node(0, label="––")
    second.add_node(1, label="a|b")
    second.add_edge(0, 1, label="m0")

    expected = nx.DiGraph()
    source1 = _matrix([[1]], cells=("shared",), muts=("a",))
    source2 = source1.copy()
    converted = object()

    trees = iter((first, second, expected))

    def fake_to_tree(value):
        tree = next(trees)
        if tree is expected:
            assert value is converted
        return tree

    monkeypatch.setattr(consensus_module.scp.ul, "to_tree", fake_to_tree)
    monkeypatch.setattr(
        consensus_module.scp.ul, "to_cfmatrix", Mock(return_value=converted)
    )

    assert consensus_module.consensus(source1, source2) is expected


def test_clt_sampler_retries_nan_softmax_and_uses_greedy_pair(monkeypatch):
    """Increase the kernel width after a numerical softmax failure."""
    calls = iter(
        [
            np.full((2, 2), np.nan),
            np.array([[0.0, 0.5], [0.5, 0.0]]),
        ]
    )
    monkeypatch.setattr(_clt_sampler, "softmax", lambda _values: next(calls))
    monkeypatch.setattr(
        _clt_sampler,
        "pairwise_distances",
        lambda _values, metric: np.array([[0.0, 1.0], [1.0, 0.0]]),
    )

    edges, probability = _clt_sampler.clt_sample_rec(
        np.array([[0.2], [0.8]]), greedy=True, c=1
    )

    assert edges == [[2, 0, 1]]
    assert probability == Decimal("0.5")


def test_clt_sampler_rejects_unimplemented_join_probability_matrix():
    """Exercise the currently unsupported cached-join-probability path."""
    with pytest.raises(AttributeError):
        _clt_sampler.clt_sample_rec(
            np.array([[0.2], [0.8]]),
            greedy=True,
            c=1,
            join_prob_matrix=np.ones((2, 2)),
        )


def test_get_samples_info_can_draw_its_own_samples(monkeypatch):
    """Return full sampling metadata when no presampled trees are supplied."""
    probability = np.array([[0.2], [0.8]])
    subtrees = [
        np.array([1, 0], dtype=np.int8),
        np.array([0, 1], dtype=np.int8),
        np.array([1, 1], dtype=np.int8),
    ]
    sampled = ([[2, 0, 1]], [subtrees], [Decimal("0.5")])
    monkeypatch.setattr(_pf, "get_samples", Mock(return_value=sampled))

    result = _pf.get_samples_info(probability, [0], 0, 1)

    assert result[2] is sampled[0]
    assert result[3] is sampled[1]
    assert result[4] is sampled[2]
    assert len(result[0]) == len(result[1]) == 1


def test_process_samples_returns_a_scalar_without_batches():
    """Return the sole weighted estimate when batches are not requested."""
    result = _pf.process_samples(
        [Decimal("0.25")],
        [Decimal("0.5")],
        [Decimal("0.25")],
    )

    assert result == Decimal("0.25")


@pytest.mark.parametrize(
    ("muts", "cells", "message"),
    [(["missing"], ["c1"], "bad muts"), (["a"], ["missing"], "bad cells")],
)
def test_partition_function_validates_requested_labels(
    muts, cells, message, capsys, monkeypatch
):
    """Reject requested mutations and cells that are absent from the matrix."""
    matrix = _matrix([[1, 0], [0, 1], [0, 0]])
    monkeypatch.setattr(scp.settings, "verbosity", "error")
    monkeypatch.setattr(scp.settings, "logfile", "")

    with pytest.raises(RuntimeError):
        _partition_function.partition_function(matrix, 0.01, 0.1, 1, 1, muts, cells)

    assert message in capsys.readouterr().out


def test_partition_function_uses_missing_value_probability(monkeypatch):
    """Convert missing genotype calls to probability one half before sampling."""
    matrix = _matrix([[3]], cells=("c1",), muts=("a",))
    seen = {}

    def fake_get_samples(probability, _n_samples):
        seen["probability"] = probability.copy()
        return [], [[np.array([1], dtype=np.int8)]], [Decimal(1)]

    monkeypatch.setattr(_partition_function, "get_samples", fake_get_samples)
    monkeypatch.setattr(
        _partition_function,
        "get_samples_info",
        Mock(return_value=([Decimal(1)], [Decimal(1)], None, None, None)),
    )
    monkeypatch.setattr(
        _partition_function, "process_samples", Mock(return_value=[Decimal(1)])
    )
    monkeypatch.setattr(
        _partition_function,
        "Parallel",
        lambda **_kwargs: lambda jobs: list(jobs),
    )
    monkeypatch.setattr(_partition_function, "delayed", lambda function: function)

    result = _partition_function.partition_function(
        matrix, 0.01, 0.1, 1, 1, ["a"], ["c1"]
    )

    assert seen["probability"].item() == 0.5
    assert result.loc["a", 0] == 1.0
