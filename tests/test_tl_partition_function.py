"""Verify partition-function estimation."""

import scphylo as scp


class TestParitionFunction:
    """Exercise mutation and cell partition-function queries."""

    def test_partition_function(self):
        """Verify that sampled partition probabilities are nonnegative."""
        df_in = scp.datasets.test()
        probs = scp.tl.partition_function(
            df_in,
            alpha=0.000001,
            beta=0.1,
            n_samples=100,
            n_batches=10,
            muts=["mut12"],
            cells=["cell6", "cell17"],
        )
        assert probs.mean(axis=1).round(4).values[0] >= 0
