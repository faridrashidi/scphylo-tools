"""Cover utility command adapters with in-process mocked work."""

from contextlib import nullcontext
from decimal import Decimal
from unittest.mock import Mock

import networkx as nx
import numpy as np
import pandas as pd

from scphylo.commands.utils import (
    _consensus,
    _partf,
    _score,
    _search,
    _trees,
    cli_utils,
)


def _invoke(command, *args):
    """Invoke a Click callback after narrowing its optional type."""
    callback = command.callback
    assert callback is not None
    return callback(*args)


def test_utils_group_preserves_command_order():
    """List utilities in their intentionally declared order."""
    assert list(cli_utils.list_commands(None)) == [
        "score",
        "search",
        "cf2newick",
        "cf2tree",
        "consensus",
        "partf",
    ]
    assert _invoke(cli_utils) is None


def test_consensus_command_reads_combines_and_writes(monkeypatch, tmp_path):
    """Delegate consensus construction and write its embedded matrix."""
    first = tmp_path / "first.CFMatrix"
    second = tmp_path / "second.CFMatrix"
    output = tmp_path / "consensus.CFMatrix"
    matrices = [object(), object()]
    final_tree = nx.DiGraph(data="consensus-matrix")
    read = Mock(side_effect=matrices)
    write = Mock()
    combine = Mock(return_value=final_tree)
    monkeypatch.setattr(_consensus.scp.io, "read", read)
    monkeypatch.setattr(_consensus.scp.io, "write", write)
    monkeypatch.setattr(_consensus.scp.tl, "consensus", combine)
    monkeypatch.setattr(
        _consensus.scp.settings, "verbosity", _consensus.scp.settings.verbosity
    )

    assert _invoke(_consensus.consensus, str(first), str(second), str(output)) is None

    combine.assert_called_once_with(*matrices)
    write.assert_called_once_with("consensus-matrix", str(output))


def test_score_command_runs_every_metric(monkeypatch):
    """Calculate and report every metric exposed by the score command."""
    ground, inferred = object(), object()
    monkeypatch.setattr(_score.scp.io, "read", Mock(side_effect=[ground, inferred]))
    monkeypatch.setattr(_score.scp.settings, "verbosity", _score.scp.settings.verbosity)
    metric_names = ("gs", "ad", "dl", "tpted", "rf", "caset", "disc", "mp3")
    metrics = {}
    for name in metric_names:
        metrics[name] = Mock(return_value=0.5)
        monkeypatch.setattr(_score.scp.tl, name, metrics[name])
    metrics["mltd"] = Mock(return_value={"normalized_similarity": 0.5})
    monkeypatch.setattr(_score.scp.tl, "mltd", metrics["mltd"])
    info = Mock()
    monkeypatch.setattr(_score.scp.logg, "info", info)

    assert _invoke(_score.score, "ground", "inferred") is None

    for metric in metrics.values():
        metric.assert_called_once_with(ground, inferred)
    assert info.call_count == 9


def test_tree_conversion_commands(monkeypatch, tmp_path):
    """Write Newick metadata and delegate clonal-tree rendering."""
    source = tmp_path / "tree.CFMatrix"
    matrix = object()
    tree = object()
    info = pd.DataFrame({"node": [1]})
    monkeypatch.setattr(_trees.scp.io, "read", Mock(return_value=matrix))
    monkeypatch.setattr(_trees.scp.ul, "to_tree", Mock(return_value=tree))
    monkeypatch.setattr(
        _trees, "_newick_info2_mutation_list", Mock(return_value=("(a);", info, []))
    )
    render = Mock()
    monkeypatch.setattr(_trees.scp.pl, "clonal_tree", render)

    assert _invoke(_trees.cf2newick, str(source)) is None
    assert _invoke(_trees.cf2tree, str(source)) is None

    assert (tmp_path / "tree.newick").read_text() == "(a);\n"
    assert (tmp_path / "tree.info2").is_file()
    render.assert_called_once_with(tree, output_file=str(tmp_path / "tree.png"))


def test_run_scistree_returns_likelihood_and_parameters(monkeypatch):
    """Return a grid-search result together with its originating rates."""
    inferred = object()
    monkeypatch.setattr(_search.scp.tl, "scistree", Mock(return_value=inferred))
    monkeypatch.setattr(_search.scp.ul, "calc_nll_matrix", Mock(return_value=1.25))

    assert _search.run_scistree("input", 0.01, 0.2, "unused") == (
        1.25,
        inferred,
        0.01,
        0.2,
    )


def test_search_command_uses_all_cores_and_selects_best_result(monkeypatch):
    """Evaluate the twenty-point grid and persist its minimum-likelihood matrix."""
    jobs_seen = {}
    matrix = object()
    best = object()

    class ImmediateParallel:
        def __init__(self, n_jobs):
            jobs_seen["n_jobs"] = n_jobs

        def __call__(self, jobs):
            values = list(jobs)
            jobs_seen["count"] = len(values)
            return values

    def fake_run(_matrix, alpha, beta, _outfile):
        return (0 if (alpha, beta) == (0.00001, 0.4) else 1, best, alpha, beta)

    monkeypatch.setattr(_search.scp.io, "read", Mock(return_value=matrix))
    monkeypatch.setattr(
        _search.scp.settings, "verbosity", _search.scp.settings.verbosity
    )
    monkeypatch.setattr(_search.scp.settings, "logfile", _search.scp.settings.logfile)
    write = Mock()
    monkeypatch.setattr(_search.scp.io, "write", write)
    monkeypatch.setattr(_search, "run_scistree", fake_run)
    monkeypatch.setattr(_search, "Parallel", ImmediateParallel)
    monkeypatch.setattr(_search, "delayed", lambda function: function)
    monkeypatch.setattr(_search, "tqdm", lambda **_kwargs: object())
    monkeypatch.setattr(_search.scp.ul, "tqdm_joblib", lambda _bar: nullcontext())
    stat = Mock()
    monkeypatch.setattr(_search.scp.ul, "stat", stat)
    monkeypatch.setattr(_search.scp.logg, "info", Mock())

    assert _invoke(_search.search, "input.SC", -1) is None

    assert jobs_seen == {"n_jobs": 20, "count": 20}
    write.assert_called_once_with(best, "input.scistree.CFMatrix")
    stat.assert_called_once_with(matrix, best, 0.00001, 0.4, 0)


def test_partf_runs_samples_in_process(monkeypatch, tmp_path):
    """Exercise the sampling closure without starting worker processes."""
    source = tmp_path / "input.SC"
    source.touch()
    matrix = pd.DataFrame([[3]], index=["c"], columns=["m"])
    seen = {}

    class ImmediateParallel:
        def __init__(self, n_jobs):
            seen["n_jobs"] = n_jobs

        def __call__(self, jobs):
            return list(jobs)

    def fake_draw(probability, greedy, c, coef):
        seen["probability"] = probability.item()
        return [], [np.array([1], dtype=np.int8)], Decimal(1)

    monkeypatch.setattr(_partf.scp.io, "read", Mock(return_value=matrix))
    monkeypatch.setattr(_partf.scp.settings, "verbosity", _partf.scp.settings.verbosity)
    monkeypatch.setattr(_partf, "Parallel", ImmediateParallel)
    monkeypatch.setattr(_partf, "delayed", lambda function: function)
    monkeypatch.setattr(_partf, "draw_sample_clt", fake_draw)
    monkeypatch.setattr(_partf.random, "choice", lambda _letters: "a")

    assert _invoke(_partf.partf, str(source), 0.01, 0.1, 1, 1) is None

    assert seen == {"n_jobs": 1, "probability": 0.5}
    samples = tmp_path / "input.partf.samples" / ("a" * 16 + ".pkl")
    assert samples.is_file()
