"""Exercise solver adapters and helpers without starting external programs."""

from contextlib import nullcontext
from types import SimpleNamespace
from unittest.mock import Mock

import networkx as nx
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from scphylo.tl.solver import (
    _bnb,
    _cardelino,
    _dendro,
    _grmt,
    _onconem,
    _phiscs,
    _sasc,
    _sbm,
    _scistree,
)
from scphylo.tl.solver.booster import _booster, _dependencies, _subsamples


def _matrix(values=((0, 1), (1, 0))):
    """Build a small labeled genotype matrix."""
    return pd.DataFrame(values, index=["c1", "c2"], columns=["a", "b"])


def test_bnb_rejects_an_unknown_bounding_strategy():
    """Validate the public branch-and-bound strategy name immediately."""
    with pytest.raises(RuntimeError):
        _bnb.bnb(_matrix(), "unknown")


def test_bnb_small_matrix_utilities():
    """Cover loop-based intersections and the abstract strategy defaults."""
    matrix = np.array([[1, 1, 0], [0, 1, 1]], dtype=np.int8)
    expected = _bnb.calculate_column_intersections(matrix, row_by_row=True)
    assert np.array_equal(
        _bnb.calculate_column_intersections(matrix, for_loop=True), expected
    )
    assert _bnb.all_None(None, None)
    assert not _bnb.all_None(None, 1)
    assert np.array_equal(_bnb.zero_or_na(np.array([0, 1, 3]), 3), [True, False, True])

    strategy = _bnb.BoundingAlgAbstract()
    with pytest.raises(RuntimeError):
        strategy.reset(matrix)
    with pytest.raises(RuntimeError):
        strategy.get_bound(sp.lil_matrix(matrix.shape))
    assert strategy.get_name() == "BoundingAlgAbstract"
    assert strategy.get_state() is None
    strategy.set_state(None)
    assert strategy.get_extra_info() == {}
    assert strategy.get_priority(1, 2, 3) == -3
    assert strategy.get_times() == {}
    assert strategy.get_init_node() is None


def test_bnb_compact_level_two_constraints_and_model():
    """Build the compact second-level constraints and weighted SAT clauses."""
    result = _bnb.make_constraints_np_matrix(
        np.array([[0, 1], [1, 0]], dtype=np.int8),
        n_levels=2,
        na_value=3,
        compact_formulation=True,
    )
    assert result.hard_constraints[1].shape == (2, 4)

    constraints = [
        np.array([[0, 0, 5, 1]], dtype=int),
        np.array([[0, 0, 6, 7]], dtype=int),
    ]
    model = _bnb.make_twosat_model_from_np(
        constraints,
        np.array([[1]], dtype=int),
        zero_vars=[1],
        na_vars=[2],
        heuristic_setting=[False, False, False, False, False],
        compact_formulation=True,
    )
    try:
        assert model.compute() is not None
    finally:
        model.delete()


@pytest.mark.parametrize(
    ("version", "expected"),
    [(1, 6), (2, 5), (3, 3), (4, 4), (5, 1), (6, 3), (7, 0)],
)
def test_twosat_priority_versions(version, expected):
    """Calculate every supported branch-and-bound queue priority."""
    strategy = _bnb.TwoSatBounding(priority_version=version)

    assert strategy.get_priority(1, 2, 3) == expected


def test_twosat_strategy_names_and_rejects_unknown_priority():
    """Encode strategy settings in its name and reject unsupported priorities."""
    strategy = _bnb.TwoSatBounding(priority_version=8)

    assert strategy.get_name().startswith("TwoSatBounding_8_")
    with pytest.raises(RuntimeError):
        strategy.get_priority(1, 2, 3)


def test_bnb_conflict_free_and_unknown_solver_paths(monkeypatch):
    """Handle conflict-free branches and an inconclusive pybnb result."""
    assert _bnb.is_conflict_free_gusfield_and_get_two_columns_in_coflicts(
        np.array([[1, 1], [1, 0]], dtype=np.int8), 3
    ) == (True, (None, None))

    problem = object.__new__(_bnb.BnB)
    problem.icf = True
    assert list(problem.branch()) == []

    result = SimpleNamespace(solution_status="unknown", termination_condition="limit")
    solver = Mock()
    solver.solve.return_value = result
    solver_factory = Mock(return_value=solver)
    monkeypatch.setattr(_bnb.pybnb.solver, "Solver", solver_factory)
    monkeypatch.setattr(_bnb, "BnB", Mock(return_value="problem"))

    output, condition = _bnb.bnb_solve(
        np.ones((2, 2), dtype=np.int8), Mock(), na_value=3, time_limit=1
    )

    assert np.array_equal(output, np.zeros((1, 1)))
    assert condition == "limit"


def test_grmt_parses_combined_labels_and_paths(monkeypatch, tmp_path):
    """Reconstruct genotypes from mocked GRMT DOT labels and paths."""
    matrix = pd.DataFrame([[0], [0]], index=["cell1", "cell2"], columns=["mutA mutB"])
    graph = nx.DiGraph()
    graph.add_node("0", label='"root"')
    graph.add_node("1", label='"mutA mutB"')
    graph.add_node("2", label='"cell1"')
    graph.add_node("3", label='"cell2"')
    graph.add_edges_from([("0", "1"), ("1", "2"), ("0", "3")])
    tempdir = SimpleNamespace(name=str(tmp_path), cleanup=Mock())
    monkeypatch.setattr(_grmt.scp.ul, "executable", Mock(return_value="grmt"))
    monkeypatch.setattr(_grmt.scp.ul, "tmpdirsys", Mock(return_value=tempdir))
    monkeypatch.setattr(_grmt.os, "system", Mock(return_value=0))
    monkeypatch.setattr(_grmt.nx.drawing.nx_pydot, "read_dot", Mock(return_value=graph))
    monkeypatch.setattr(_grmt.scp.ul, "stat", Mock())

    result = _grmt.grmt(matrix, 0.01, 0.2, n_iters=1, n_threads=1)

    assert result.loc["cell1", "mutA mutB"] == 1
    assert result.loc["cell2", "mutA mutB"] == 0
    tempdir.cleanup.assert_called_once_with()


def test_sasc_placeholder_returns_none():
    """Document the current SASC placeholder behavior."""
    assert _sasc.sasc() is None


@pytest.mark.parametrize(
    ("module", "importer", "function", "arguments"),
    [
        (_cardelino, "import_rpy2", _cardelino.cardelino, (object(), "free")),
        (_dendro, "import_rpy2", _dendro.dendro, (object(),)),
        (_onconem, "import_rpy2", _onconem.onconem, (_matrix(), 0.01, 0.2)),
        (_sbm, "import_graph_tool", _sbm.sbm, (_matrix(),)),
    ],
)
def test_optional_solvers_report_missing_packages(
    monkeypatch, module, importer, function, arguments
):
    """Fail clearly when an optional R or graph-tool package is unavailable."""
    monkeypatch.setattr(module.scp.ul, importer, Mock(return_value=(None, True)))

    with pytest.raises(RuntimeError):
        function(*arguments)


@pytest.mark.parametrize(
    ("function", "arguments"),
    [
        (_phiscs.phiscsi, (_matrix(), 0.01, 0.2)),
        (_phiscs.phiscsi_bulk, (_matrix(), 0.01, 0.2)),
        (_phiscs.phiscs_readcount, (object(), 0.01, 0.2)),
    ],
)
def test_gurobi_solvers_report_a_missing_package(monkeypatch, function, arguments):
    """Fail before model construction when Gurobi is unavailable."""
    monkeypatch.setattr(
        _phiscs.scp.ul, "import_gurobi", Mock(return_value=(None, True))
    )

    with pytest.raises(RuntimeError):
        function(*arguments)


def test_phiscs_binary_zero_alpha_and_bulk_modes(monkeypatch):
    """Cover one-sided correction and both mutation-elimination modes."""
    monkeypatch.setattr(_phiscs.scp.ul, "stat", Mock())
    mixed = _matrix()
    no_zeros = _matrix(((1, 1), (1, 1)))

    outputs = [
        _phiscs.phiscsb(mixed, 0, 0.2, experiment=True),
        _phiscs.phiscsb_bulk(no_zeros, 0, 0.2, kmax=1),
        _phiscs.phiscsb_bulk(mixed, 0.1, 0.2, kmax=0),
    ]

    assert all(output.shape == mixed.shape for output in outputs)


class _ReadCountData:
    """Minimal read-count input used to reach rScisTree mode validation."""

    obs_names = pd.Index(["c"])
    var_names = pd.Index(["m"])
    layers = {"mutant": np.array([[1]]), "total": np.array([[2]])}

    @staticmethod
    def to_df():
        """Return a one-call genotype matrix."""
        return pd.DataFrame([[1]], index=["c"], columns=["m"])


def test_rscistree_ternary_and_invalid_modes(monkeypatch, tmp_path):
    """Select ternary ScisTree mode and reject an invalid mode."""
    tempdir = SimpleNamespace(name=str(tmp_path), cleanup=Mock())
    monkeypatch.setattr(_scistree.scp.ul, "tmpdirsys", Mock(return_value=tempdir))

    class StopAfterModeSelection(Exception):
        pass

    run_scprob = Mock(side_effect=StopAfterModeSelection)
    monkeypatch.setattr(_scistree, "run_scprob", run_scprob)
    with pytest.raises(StopAfterModeSelection):
        _scistree.rscistree(_ReadCountData(), mode="ternary")
    assert run_scprob.call_args.args[0][-1] == "1"

    with pytest.raises(RuntimeError):
        _scistree.rscistree(_ReadCountData(), mode="invalid")


def test_iscistree_handles_tiny_distance_matrices(monkeypatch):
    """Cover neighbor joining's one- and two-taxon special cases."""
    monkeypatch.setattr(_scistree.scp.ul, "stat", Mock())
    one_cell = pd.DataFrame([[1]], index=["c"], columns=["m"])
    result = _scistree.iscistree(one_cell, 0.01, 0.2, n_iters=1)
    assert result.loc["c", "m"] == 1

    no_cells = pd.DataFrame(np.empty((0, 1), dtype=int), columns=["m"])
    with pytest.raises(ValueError):
        _scistree.iscistree(no_cells, 0.01, 0.2, n_iters=1)


def test_booster_noop_modes_and_sample_size_validation(monkeypatch, tmp_path):
    """Return no reconstruction and expose the invalid automatic sample-size path."""
    matrix = _matrix()
    monkeypatch.setattr(_booster.scp.ul, "mkdir", Mock(return_value=str(tmp_path)))

    assert (
        _booster.booster(
            matrix.copy(),
            0.01,
            0.2,
            subsample_dir=str(tmp_path),
            no_subsampling=True,
            no_dependencies=True,
            no_reconstruction=True,
        )
        is None
    )
    with pytest.raises(TypeError):
        _booster.booster(
            matrix.copy(),
            0.01,
            0.2,
            sample_size=None,
            subsample_dir=str(tmp_path),
            no_subsampling=True,
            no_dependencies=True,
            no_reconstruction=True,
        )


@pytest.mark.parametrize(
    ("sample_on", "solver", "expected_error"),
    [("invalid", "scite", RuntimeError), ("muts", "invalid", RuntimeError)],
)
def test_subsampling_validates_modes(
    monkeypatch, tmp_path, sample_on, solver, expected_error
):
    """Reject unknown sampling dimensions and solver names in process."""
    monkeypatch.setattr(_subsamples.scp.ul, "with_timeout", lambda _limit: lambda f: f)
    monkeypatch.setattr(_subsamples, "delayed", lambda function: function)
    monkeypatch.setattr(
        _subsamples, "Parallel", lambda **_kwargs: lambda jobs: list(jobs)
    )
    monkeypatch.setattr(_subsamples, "tqdm", lambda **_kwargs: object())
    monkeypatch.setattr(_subsamples.scp.ul, "tqdm_joblib", lambda _bar: nullcontext())

    with pytest.raises(expected_error):
        _subsamples.subsampling(
            _matrix(),
            0.01,
            0.2,
            solver,
            sample_on,
            1,
            1,
            0,
            1,
            1,
            1,
            str(tmp_path),
            True,
        )


def test_subsampling_skips_too_small_scistree_input(monkeypatch, tmp_path):
    """Skip ScisTree when filtering leaves fewer than two mutations."""
    monkeypatch.setattr(_subsamples.scp.ul, "with_timeout", lambda _limit: lambda f: f)
    monkeypatch.setattr(_subsamples, "delayed", lambda function: function)
    monkeypatch.setattr(
        _subsamples, "Parallel", lambda **_kwargs: lambda jobs: list(jobs)
    )
    monkeypatch.setattr(_subsamples, "tqdm", lambda **_kwargs: object())
    monkeypatch.setattr(_subsamples.scp.ul, "tqdm_joblib", lambda _bar: nullcontext())
    solve = Mock()
    monkeypatch.setattr(_subsamples.scp.tl, "scistree", solve)

    assert (
        _subsamples.subsampling(
            _matrix(((0, 0), (0, 0))),
            0.01,
            0.2,
            "scistree",
            "muts",
            1,
            1,
            0,
            1,
            1,
            1,
            str(tmp_path),
            True,
        )
        is None
    )
    solve.assert_not_called()


@pytest.mark.parametrize(
    ("values", "expected"),
    [
        ({"c": {"a": 0, "b": 1}}, _dependencies.UNDEFINED_DEPENDENCY),
        ({"c": {"a": 1, "b": 0}}, _dependencies.UNDEFINED_DEPENDENCY),
        ({"c": {"a": 0, "b": 0}}, _dependencies.UNDEFINED_DEPENDENCY),
        ({"c": {"a": 2, "b": 2}}, _dependencies.UNDEFINED_DEPENDENCY),
    ],
)
def test_dependency_undefined_cases(values, expected):
    """Classify one-sided, empty, and ignored nonbinary observations."""
    assert (
        _dependencies.get_dependency_from_conflict_free_matrix("a", "b", values)
        == expected
    )


def test_dependency_number_format_helpers():
    """Format numeric values and leave nonnumeric values untouched."""
    assert _dependencies.is_float("1.2")
    assert not _dependencies.is_float("bad")
    assert _dependencies.float_to_string(1.2) == "1.20"
    assert _dependencies.float_to_string("bad") == "bad"
