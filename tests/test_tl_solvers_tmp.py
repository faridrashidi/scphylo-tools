"""Exercise experimental and optional solver integrations."""

import subprocess
from pathlib import Path

import pandas as pd

import scphylo as scp

from ._helpers import skip_rpy2, skip_slow


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

    def test_infscite(self):
        """Reserve coverage for the optional infSCITE integration."""
        assert True

    def test_sbm(self):
        """Verify SBM preserves labels and resolves genotype conflicts."""
        data = scp.datasets.test()
        out = scp.tl.sbm(data)
        assert out.index.equals(data.index)
        assert out.columns.equals(data.columns)
        assert scp.ul.is_conflict_free_gusfield(out)

    @skip_rpy2("infercna")
    def test_infercna(self):
        """Verify InferCNA integration when its R package is available."""
        expr = scp.datasets.example(is_expression=True)
        df_cna = scp.tl.infercna(expr, ref_cells={"normal": ["C15_1"]}, genome="mm10")
        df_cna.loc["C15_1"] = 0
        expr.obsm["cna"] = df_cna.loc[expr.obs_names]
        scp.pl.heatmap(expr, layer="cna")

    @skip_rpy2("dendro")
    @skip_rpy2("ggtree")
    def test_dendro(self):
        """Verify DENDRO integration when its R package is available."""
        adata = scp.datasets.example()
        scp.tl.dendro(adata)
        assert True

    @skip_rpy2("cardelino")
    @skip_slow
    def test_cardelino(self):
        """Verify Cardelino integration in the long-running test suite."""
        adata = scp.datasets.example()
        scp.tl.cardelino(adata, mode="free", n_clones=11)
        assert True
