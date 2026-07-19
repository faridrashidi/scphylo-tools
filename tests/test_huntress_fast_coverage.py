"""Fast unit coverage for the Huntress reconstruction internals."""

from importlib import import_module

import numpy as np
import pandas as pd

huntress = import_module("scphylo.tl.solver.huntress._huntress")


class _Queue:
    """Provide the small queue interface used by Huntress."""

    def __init__(self):
        self.items = []

    def put(self, item):
        """Append an item."""
        self.items.append(item)

    def get(self):
        """Return the oldest item."""
        return self.items.pop(0)


def test_fixed_tuning_runs_greedy_na_without_processes(monkeypatch):
    """Exercise fixed-parameter reconstruction and both greedy assignments."""
    matrix = np.array([[1, 1, 0], [1, 0, 1], [0, 1, 1]])
    frame = pd.DataFrame(matrix)
    monkeypatch.setattr(huntress, "Queue", _Queue)
    monkeypatch.setattr(
        huntress,
        "postprocess_col",
        lambda _frame, reconstructed, **_kwargs: reconstructed,
    )

    reconstructed, running_time = huntress.Reconstruct(
        frame,
        auto_tune=0,
        overlapp_coeff=0.8,
        hist_coeff=20,
    )

    np.testing.assert_array_equal(
        reconstructed,
        [[1, 1, 0], [0, 0, 1], [0, 0, 1]],
    )
    assert running_time >= 0


def test_auto_tuning_chunks_and_selects_later_result(monkeypatch):
    """Cover uneven tuning chunks and selection of a later lower score."""

    class _Process:
        def __init__(self, target, args):
            self.target = target
            self.args = args

        def start(self):
            self.target(*self.args)

        def join(self):
            return None

    def auto_tune(queue, tuning_range, *_args):
        score = len(tuning_range)
        queue.put([np.full((2, 2), score == 66), score])

    monkeypatch.setattr(huntress, "Queue", _Queue)
    monkeypatch.setattr(huntress, "Process", _Process)
    monkeypatch.setattr(huntress, "Auto_fnfp", auto_tune)
    monkeypatch.setattr(
        huntress,
        "postprocess_col",
        lambda _frame, reconstructed, **_kwargs: reconstructed,
    )

    reconstructed, running_time = huntress.Reconstruct(
        pd.DataFrame(np.zeros((2, 2), dtype=int)),
        n_proc=3,
    )

    np.testing.assert_array_equal(reconstructed, np.ones((2, 2), dtype=bool))
    assert running_time >= 0


def test_auto_fnfp_handles_empty_range_and_keeps_best_candidate(monkeypatch):
    """Cover the tuning worker's empty guard and improving candidates."""
    matrix_input = np.array([[1, 0], [0, 1]], dtype=bool)
    matrix_raw = np.array([[1, 3], [0, 1]])
    estimated = matrix_input.astype(float)
    queue = _Queue()

    huntress.Auto_fnfp(queue, [], matrix_input, estimated, matrix_raw, 10, 1)
    assert queue.items == []

    candidates = iter(
        [
            np.zeros((2, 2), dtype=bool),
            np.ones((2, 2), dtype=bool),
            matrix_input,
        ]
    )
    monkeypatch.setattr(
        huntress,
        "greedyPtreeNA",
        lambda *_args: [None, None, next(candidates)],
    )

    huntress.Auto_fnfp(
        queue,
        [(1, 0.1), (3, 0.2)],
        matrix_input,
        estimated,
        matrix_raw,
        10,
        1,
    )

    reconstructed, distance = queue.get()
    np.testing.assert_array_equal(reconstructed, matrix_input)
    assert distance == 0


def test_conversion_helpers_cover_missing_only_column():
    """Convert raw values and estimate a mutation containing only missing data."""
    frame = pd.DataFrame([[3, 0], [3, 1]])

    np.testing.assert_array_equal(
        huntress.ReadFfile(frame),
        [[True, False], [True, True]],
    )
    np.testing.assert_array_equal(
        huntress.Estimated_Matrix(frame),
        [[0.0, 0.0], [0.0, 1.0]],
    )


def test_postprocess_retains_an_improving_iteration(monkeypatch):
    """Retain one improved column assignment, then stop when it stabilizes."""
    assignments = iter(
        [
            np.zeros((2, 2), dtype=bool),
            np.ones((2, 2), dtype=bool),
            np.ones((2, 2), dtype=bool),
        ]
    )
    monkeypatch.setattr(
        huntress,
        "c_m_col",
        lambda *_args, **_kwargs: next(assignments),
    )
    monkeypatch.setattr(
        huntress,
        "c_m_row",
        lambda _input, nodes, **_kwargs: nodes,
    )

    reconstructed = huntress.postprocess_col(
        pd.DataFrame(np.ones((2, 2), dtype=int)),
        np.zeros((2, 2), dtype=bool),
        pfn=0.1,
        pfp=0.01,
    )

    np.testing.assert_array_equal(reconstructed, np.ones((2, 2), dtype=bool))
