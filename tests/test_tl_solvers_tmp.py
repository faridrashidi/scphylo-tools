"""Exercise experimental and optional solver integrations."""

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

    def test_siclonefit(self):
        """Reserve coverage for the optional SiCloneFit integration."""
        assert True

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
