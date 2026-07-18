import pytest

import scphylo as scp


class TestSolversTmp:
    def test_titch_algorithm(self):
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
        adata = scp.datasets.colorectal2(readcount=True)
        df_out = scp.tl.rscistree(adata, mode="haploid")
        is_cf = scp.ul.is_conflict_free_gusfield(df_out)
        assert is_cf

    def test_iscistree(self):
        df_in = scp.datasets.test()
        df_out = scp.tl.iscistree(df_in, alpha=0.0000001, beta=0.1)
        is_cf = scp.ul.is_conflict_free_gusfield(df_out)
        assert is_cf

    def test_siclonefit(self):
        assert True

    def test_infscite(self):
        assert True

    def test_sbm(self):
        data = scp.datasets.test()
        out = scp.tl.sbm(data)
        assert out.index.equals(data.index)
        assert out.columns.equals(data.columns)
        assert scp.ul.is_conflict_free_gusfield(out)

    @pytest.mark.skip(reason="Unable to import a package on GitHub!")
    def test_infercna(self):
        expr = scp.datasets.example(is_expression=True)
        df_cna = scp.tl.infercna(expr, ref_cells={"normal": ["C15_1"]}, genome="mm10")
        df_cna.loc["C15_1"] = 0
        expr.obsm["cna"] = df_cna.loc[expr.obs_names]
        scp.pl.heatmap(expr, layer="cna")

    @pytest.mark.skip(reason="Unable to import a package on GitHub!")
    def test_dendro(self):
        adata = scp.datasets.example()
        scp.tl.dendro(adata)
        assert True

    @pytest.mark.skip(reason="Takes 6 minutes!")
    def test_cardelino(self):
        adata = scp.datasets.example()
        scp.tl.cardelino(adata, mode="free", n_clones=11)
        assert True
