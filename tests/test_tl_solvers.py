"""Verify phylogeny solver implementations and booster backends."""

import builtins
from types import SimpleNamespace

import numpy as np
from anndata import AnnData

import scphylo as scp

from ._helpers import fake_gurobi


def _use_fake_gurobi(monkeypatch, value_for_variable=None):
    """Install the deterministic solver double used by PhISCS adapter tests."""
    if value_for_variable is None:

        def value_for_variable(model, index, name):
            return int(name is not None and name.startswith("Y(0,"))

    gurobi = fake_gurobi(value_for_variable)
    monkeypatch.setattr(scp.ul, "import_gurobi", lambda: (gurobi, False))
    return gurobi


def _use_fake_onconem(monkeypatch, shape):
    """Install a fast OncoNEM API double while retaining real R conversion."""
    import rpy2.robjects as ro
    import rpy2.robjects.packages as rpackages

    class OncoNEM:
        def __init__(self, data_shape):
            self.shape = data_shape

        def oncoNEM(self, Data, **kwargs):
            return ro.r("list(search=function(delta) NULL)")

        def expandOncoNEM(self, onem, **kwargs):
            return onem

        def clusterOncoNEM(self, **kwargs):
            return ro.ListVector({"g": ro.IntVector([0]), "clones": ro.IntVector([0])})

        def oncoNEMposteriors(self, **kwargs):
            n_muts, n_cells = self.shape
            values = np.zeros((n_muts, n_cells), dtype=int)
            values[:, 0] = 1
            p_mut = ro.r.matrix(
                ro.IntVector(values.ravel(order="F")),
                nrow=n_muts,
                ncol=n_cells,
            )
            return ro.ListVector({"p_mut": p_mut})

    backend = OncoNEM(shape)
    monkeypatch.setattr(scp.ul, "import_rpy2", lambda *args, **kwargs: (backend, False))
    monkeypatch.setattr(
        rpackages,
        "importr",
        lambda name: SimpleNamespace(get_edgelist=lambda graph: None),
    )


class TestSolvers:
    """Exercise solver APIs on a shared genotype matrix."""

    def setup_method(self):
        """Load a fresh genotype matrix for each solver test."""
        self.df_in = scp.datasets.test()

    def test_scite(self):
        """Verify that SCITE produces a conflict-free matrix."""
        df_out = scp.tl.scite(
            self.df_in, alpha=0.0000001, beta=0.1, n_restarts=3, n_iters=1000
        )
        assert scp.ul.is_conflict_free_gusfield(df_out)
        assert scp.ul.is_conflict_free(df_out)

    def test_bnb_simulated(self):
        """Verify branch-and-bound with simulated-data bounds."""
        df_out = scp.tl.bnb(self.df_in, bounding="simulated")
        assert scp.ul.is_conflict_free_gusfield(df_out)

    def test_bnb_serial_does_not_import_mpi(self, monkeypatch):
        """Verify that serial branch-and-bound does not require mpi4py."""
        original_import = builtins.__import__

        def reject_mpi_import(name, *args, **kwargs):
            if name == "mpi4py" or name.startswith("mpi4py."):
                raise AssertionError("serial BnB attempted to import mpi4py")
            return original_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", reject_mpi_import)
        df_out = scp.tl.bnb(
            self.df_in,
            bounding="simulated",
            time_limit=20,
        )
        assert scp.ul.is_conflict_free_gusfield(df_out)

    def test_bnb_real(self):
        """Verify branch-and-bound with real-data bounds."""
        df_out = scp.tl.bnb(self.df_in, bounding="real", time_limit=20)
        assert scp.ul.is_conflict_free_gusfield(df_out)

    def test_scistree(self):
        """Verify that ScisTree produces a conflict-free matrix."""
        df_out = scp.tl.scistree(self.df_in, alpha=0.0000001, beta=0.1)
        assert scp.ul.is_conflict_free_gusfield(df_out)

    def test_phiscsb(self):
        """Verify that binary PhISCS produces a conflict-free matrix."""
        df_out = scp.tl.phiscsb(self.df_in, alpha=0.0000001, beta=0.1)
        assert scp.ul.is_conflict_free_gusfield(df_out)

    def test_huntress_both(self):
        """Verify HUNTRESS correction with both error types enabled."""
        df_out = scp.tl.huntress(self.df_in, alpha=0.0000001, beta=0.1)
        assert scp.ul.is_conflict_free_gusfield(df_out)

    def test_huntress_fn(self):
        """Verify HUNTRESS correction with false negatives only."""
        df_out = scp.tl.huntress(self.df_in, alpha=0, beta=0.1)
        assert scp.ul.is_conflict_free_gusfield(df_out)

    def test_onconem(self, monkeypatch):
        """Verify OncoNEM conversion, search plumbing, and output labels."""
        df_in = self.df_in.iloc[:2, :3]
        _use_fake_onconem(monkeypatch, (df_in.shape[1], df_in.shape[0]))
        df_out = scp.tl.onconem(df_in, alpha=0.0000001, beta=0.1)
        assert scp.ul.is_conflict_free_gusfield(df_out)
        assert df_out.index.tolist() == df_in.index.tolist()
        assert df_out.columns.tolist() == df_in.columns.tolist()

    def test_phiscsi(self, monkeypatch):
        """Verify integer PhISCS model construction and output labels."""
        gurobi = _use_fake_gurobi(monkeypatch)
        df_in = self.df_in.iloc[:2, :3]
        df_out = scp.tl.phiscsi(df_in, alpha=0.0000001, beta=0.1)
        assert scp.ul.is_conflict_free_gusfield(df_out)
        assert df_out.index.tolist() == df_in.index.tolist()
        assert df_out.columns.tolist() == df_in.columns.tolist()
        assert gurobi.models[0].Params.Threads == 1

    def test_phiscsi_bulk_1(self, monkeypatch):
        """Verify bulk-aware PhISCS model construction with VAF information."""
        gurobi = _use_fake_gurobi(monkeypatch)
        df_in = self.df_in.iloc[:2, :2]
        vaf_info = df_in.T.iloc[:, :1].astype(float)
        vaf_info.iloc[:, 0] = [0.2, 0.4]
        df_out = scp.tl.phiscsi_bulk(
            df_in,
            alpha=0.01,
            beta=0.1,
            delta=0.2,
            kmax=1,
            vaf_info=vaf_info,
        )
        assert scp.ul.is_conflict_free_gusfield(df_out)
        assert df_out.shape == df_in.shape
        assert gurobi.models[0].objective[1] == gurobi.GRB.MAXIMIZE

    def test_phiscsi_bulk_2(self, monkeypatch):
        """Verify bulk-aware PhISCS without VAF information."""
        _use_fake_gurobi(monkeypatch)
        df_in = self.df_in.iloc[:2, :2]
        df_out = scp.tl.phiscsi_bulk(df_in, 0.01, 0.1, time_limit=1)
        assert scp.ul.is_conflict_free_gusfield(df_out)
        assert df_out.iloc[0].eq(1).all()
        assert df_out.iloc[1].eq(0).all()

    def test_phiscsi_bulk_3(self, monkeypatch):
        """Verify that bulk-aware PhISCS zeros eliminated mutations."""

        def solution_values(model, index, name):
            return int(name is not None and (name.startswith("Y(") or name == "K[1]"))

        _use_fake_gurobi(monkeypatch, solution_values)
        df_in = self.df_in.iloc[:2, :2]
        df_out = scp.tl.phiscsi_bulk(df_in, 0.01, 0.1, kmax=1, time_limit=1)
        assert df_out.iloc[:, 0].eq(1).all()
        assert df_out.iloc[:, 1].eq(0).all()

    def test_phiscs_readcount(self, monkeypatch):
        """Verify fast PhISCS read-count model construction."""
        _use_fake_gurobi(monkeypatch)
        adata = AnnData(
            np.array([[0, 1], [1, 0]], dtype=np.int8),
            layers={
                "mutant": np.array([[0, 6], [4, 0]], dtype=np.int8),
                "total": np.full((2, 2), 10, dtype=np.int8),
            },
        )
        adata.obs_names = ["cell-a", "cell-b"]
        adata.var_names = ["mut-a", "mut-b"]
        df_out = scp.tl.phiscs_readcount(adata, alpha=0.01, beta=0.19)
        assert scp.ul.is_conflict_free_gusfield(df_out)
        assert df_out.index.tolist() == adata.obs_names.tolist()
        assert df_out.columns.tolist() == adata.var_names.tolist()

    def test_booster_phiscs(self):
        """Verify booster reconstruction with the PhISCS backend."""
        df_out = scp.tl.booster(
            self.df_in,
            alpha=0.0000001,
            beta=0.1,
            solver="PhISCS",
            sample_on="muts",
            sample_size=10,
            n_samples=20,
            begin_index=0,
            n_jobs=1,
            time_limit=120,
            dep_weight=5,
        )
        assert scp.ul.is_conflict_free_gusfield(df_out)

    def test_booster_scite(self):
        """Verify booster reconstruction with the SCITE backend."""
        df_out = scp.tl.booster(
            self.df_in,
            alpha=0.0000001,
            beta=0.1,
            solver="SCITE",
            sample_on="muts",
            sample_size=10,
            n_samples=20,
            begin_index=0,
            n_jobs=1,
            n_iterations=10000,
            dep_weight=5,
        )
        assert scp.ul.is_conflict_free_gusfield(df_out)

    def test_booster_scistree_on_cells(self):
        """Verify cell-sampled booster reconstruction with ScisTree."""
        df_out = scp.tl.booster(
            self.df_in,
            alpha=0.0000001,
            beta=0.1,
            solver="ScisTree",
            sample_on="cells",
            sample_size=10,
            n_samples=20,
            begin_index=0,
            n_jobs=1,
            n_iterations=10000,
            dep_weight=5,
        )
        assert scp.ul.is_conflict_free_gusfield(df_out)
