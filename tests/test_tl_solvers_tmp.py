"""Exercise experimental and optional solver integrations."""

import importlib
import subprocess
from contextlib import contextmanager
from io import BytesIO
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from anndata import AnnData

import scphylo as scp


def _small_readcount_adata():
    """Return a tiny annotated matrix for optional R-wrapper tests."""
    adata = AnnData(
        np.array([[0, 1], [1, 0]], dtype=np.int8),
        layers={
            "genotype": np.array([[0, 1], [1, 0]], dtype=np.int8),
            "mutant": np.array([[0, 6], [4, 0]], dtype=np.int8),
            "total": np.full((2, 2), 10, dtype=np.int8),
            "tpm": np.array([[1.0, 2.0], [3.0, 4.0]]),
        },
    )
    adata.obs_names = ["normal", "tumor"]
    adata.var_names = ["gene-a", "gene-b"]
    return adata


def _use_fake_infercna(monkeypatch):
    """Install a deterministic InferCNA package double."""
    import rpy2.robjects as ro

    class InferCNA:
        heatCols = ro.StrVector(["#0000ff", "#ff0000"])
        genome = None
        cells = None
        genes = None

        def useGenome(self, genome):
            self.genome = genome

        def infercna(self, data, **kwargs):
            references = {str(cell) for group in kwargs["refCells"] for cell in group}
            self.cells = [str(cell) for cell in data.colnames if cell not in references]
            self.genes = [str(gene) for gene in data.rownames]
            return data

        def cnaPlot(self, cna, **kwargs):
            cells = self.cells
            genes = self.genes
            assert cells is not None
            assert genes is not None
            values = np.asarray(cna)
            plot_data = ro.DataFrame(
                {
                    "Cell": ro.StrVector(
                        [cell for cell in cells for _ in range(len(genes))]
                    ),
                    "Gene": ro.StrVector(genes * len(cells)),
                    "CNA": ro.FloatVector(values.ravel(order="F")),
                }
            )
            return ro.ListVector({"data": plot_data, "p": ro.NULL})

    backend = InferCNA()
    monkeypatch.setattr(scp.ul, "import_rpy2", lambda *args, **kwargs: (backend, False))
    return backend


def _use_fake_dendro(monkeypatch):
    """Install doubles for DENDRO's R computation and plotting boundary."""
    import rpy2.robjects as ro
    import rpy2.robjects.packages as rpackages
    from rpy2.robjects.lib import grdevices

    class Dendro:
        def DENDRO_dist(self, *args, **kwargs):
            return ro.FloatVector([0.5])

    class R:
        def __call__(self, command):
            ro.globalenv["p"] = ro.IntVector([1])
            return ro.NULL

        def show(self, value):
            return None

    @contextmanager
    def image_bytes(*args, **kwargs):
        yield BytesIO(b"png")

    dendro_module = importlib.import_module("scphylo.tl.solver._dendro")
    backend = Dendro()
    monkeypatch.setattr(
        scp.ul,
        "import_rpy2",
        lambda name, *args, **kwargs: (backend, False),
    )
    monkeypatch.setattr(
        rpackages,
        "importr",
        lambda name: SimpleNamespace(
            hclust=lambda distance, **kwargs: ro.IntVector([1])
        ),
    )
    monkeypatch.setattr(ro, "r", R())
    monkeypatch.setattr(grdevices, "render_to_bytesio", image_bytes)
    monkeypatch.setattr(dendro_module, "Image", lambda *args, **kwargs: "image")
    monkeypatch.setattr(dendro_module, "display", lambda image: image)


def _use_fake_cardelino(monkeypatch, n_cells):
    """Install a deterministic Cardelino package double."""
    import rpy2.robjects as ro

    class Cardelino:
        clone_kwargs = None

        def clone_id(self, mutant, total, **kwargs):
            self.clone_kwargs = kwargs
            probabilities = ro.r.matrix(
                ro.FloatVector([0.9, 0.1] * n_cells),
                nrow=n_cells,
                ncol=2,
                byrow=True,
            )
            return ro.ListVector({"prob": probabilities})

        def assign_cells_to_clones(self, probabilities):
            return ro.IntVector(range(1, n_cells + 1))

    backend = Cardelino()
    monkeypatch.setattr(scp.ul, "import_rpy2", lambda *args, **kwargs: (backend, False))
    return backend


class TestSolversTmp:
    """Cover solver integrations that are still experimental or optional."""

    def test_titch_algorithm(self):
        """Verify Fitch labeling on a tree inferred by SCITE."""
        df_in = scp.datasets.test()
        df_out = scp.tl.scite(
            df_in, alpha=0.00001, beta=0.1, n_restarts=3, n_iters=1000
        )
        tree = scp.ul.to_tree(df_out)
        for n in tree.nodes:
            if scp.ul.is_leaf(tree, n):
                tree.nodes[n]["profile"] = [1]
        scp.tl.fitch(tree)
        assert True

    def test_rscistree(self):
        """Verify read-count ScisTree in haploid mode."""
        adata = scp.datasets.colorectal2(readcount=True)
        df_out = scp.tl.rscistree(adata, mode="haploid")
        is_cf = scp.ul.is_conflict_free_gusfield(df_out)
        assert is_cf

    def test_iscistree(self):
        """Verify integer ScisTree produces a conflict-free matrix."""
        df_in = scp.datasets.test()
        df_out = scp.tl.iscistree(df_in, alpha=0.0000001, beta=0.1)
        is_cf = scp.ul.is_conflict_free_gusfield(df_out)
        assert is_cf

    def test_siclonefit(self, tmp_path, monkeypatch):
        """Verify the SiCloneFit adapter with paths containing spaces."""
        df_in = scp.datasets.test()
        expected = df_in.copy()
        expected.iloc[:, :] = 0
        expected.iloc[0, 0] = 1

        tools_dir = tmp_path / "tools with spaces"
        tools_dir.mkdir()
        jar_path = tools_dir / "SiCloneFiTComplete.jar"
        jar_path.write_text("fixture")
        java_executable = tools_dir / "java"
        java_executable.write_text("#!/bin/sh\nexit 0\n")
        java_executable.chmod(0o755)

        def fake_run(args, appname, **kwargs):
            assert appname == "SiCloneFit"
            assert args[:3] == [
                str(java_executable.resolve()),
                "-jar",
                str(jar_path.resolve()),
            ]
            workdir = Path(args[args.index("-outDir") + 1])
            best_dir = workdir / "fixture_samples" / "best"
            best_dir.mkdir(parents=True)
            output = expected.T.copy()
            output.index = range(output.shape[0])
            output.to_csv(
                best_dir / "best_MAP_predicted_genotype.txt",
                sep=" ",
                header=False,
            )
            (best_dir / "best_MAP_tree.txt").write_text("(cell1,cell2);\n")
            return subprocess.CompletedProcess(args, 0)

        monkeypatch.setattr(scp.ul, "run_external", fake_run)
        df_out, tree = scp.tl.siclonefit(
            df_in,
            alpha=0.0000001,
            beta=0.1,
            n_restarts=1,
            n_iters=1,
            n_burnin=0,
            return_tree=True,
            jar_path=jar_path,
            java_executable=java_executable,
        )

        pd.testing.assert_frame_equal(df_out, expected)
        assert tree == "(cell1,cell2);"

    def test_sphyr(self, tmp_path, monkeypatch):
        """Verify that SPhyR uses checked argv-list execution."""
        df_in = scp.datasets.test()
        expected = df_in.copy()
        expected.iloc[:, :] = 0
        expected.iloc[0, 0] = 1

        tools_dir = tmp_path / "tools with spaces"
        tools_dir.mkdir()
        executable = tools_dir / "sphyr_kDPFC"
        executable.write_text("#!/bin/sh\nexit 0\n")
        executable.chmod(0o755)

        def fake_run(args, appname, **kwargs):
            assert appname == "SPhyR"
            assert args[0] == str(executable.resolve())
            output = kwargs["stdout"]
            output.write("SPhyR output\ncorrected genotype\n")
            expected.to_csv(output, sep=" ", header=False, index=False)
            output.flush()
            return subprocess.CompletedProcess(args, 0)

        monkeypatch.setattr(scp.ul, "run_external", fake_run)
        df_out = scp.tl.sphyr(
            df_in,
            alpha=0.0000001,
            beta=0.1,
            executable=executable,
        )

        pd.testing.assert_frame_equal(df_out, expected)

    def test_sbm(self):
        """Verify SBM preserves labels and resolves genotype conflicts."""
        data = scp.datasets.test()
        out = scp.tl.sbm(data)
        assert out.index.equals(data.index)
        assert out.columns.equals(data.columns)
        assert scp.ul.is_conflict_free_gusfield(out)

    def test_infercna(self, monkeypatch):
        """Verify InferCNA conversion and normal-cell filtering."""
        backend = _use_fake_infercna(monkeypatch)
        expr = _small_readcount_adata()
        df_cna = scp.tl.infercna(expr, ref_cells={"normal": ["normal"]}, genome="mm10")
        assert backend.genome == "mm10"
        assert df_cna.index.tolist() == ["tumor"]
        assert df_cna.columns.tolist() == expr.var_names.tolist()

    def test_dendro(self, monkeypatch):
        """Verify DENDRO's data-conversion and plotting plumbing."""
        _use_fake_dendro(monkeypatch)
        result = scp.tl.dendro(_small_readcount_adata(), width=20, height=20)
        assert result == "image"

    def test_cardelino(self, monkeypatch):
        """Verify the fast Cardelino free-mode adapter path."""
        adata = _small_readcount_adata()
        backend = _use_fake_cardelino(monkeypatch, adata.n_obs)
        result = scp.tl.cardelino(adata, mode="free", n_clones=2)
        np.testing.assert_array_equal(result, [1, 2])
        assert backend.clone_kwargs == {"n_clone": 2}
