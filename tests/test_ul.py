import numpy as np
import pandas as pd
import pytest

import scphylo as scp

from ._helpers import skip_graphviz


class TestUtils:
    def test_hclustering_1(self):
        np.random.seed(0)
        df = pd.DataFrame(np.random.randint(0, 2, size=(10, 10)))
        clusters = scp.ul.hclustering(df)
        assert clusters[6].value_counts()[5] == 3

    def test_hclustering_2(self):
        adata = scp.datasets.example()
        clusters = scp.ul.hclustering(adata.to_df(), metric="l1")
        assert len(clusters) == 81

    def test_hclustering_3(self):
        adata = scp.datasets.example()
        clusters = scp.ul.hclustering(adata.to_df(), metric="cosine")
        assert len(clusters) == 81

    def test_dist_dendro(self):
        scp.settings.verbosity = 0
        adata = scp.datasets.example()
        dist = scp.ul.dist_dendro(adata)
        assert dist.sum() > 262456
        assert dist.sum() < 262457

    def test_tree_to_cfmatrix(self):
        df_in = scp.datasets.test()
        df_out = scp.tl.phiscsb(df_in, alpha=0.0000001, beta=0.1)
        tree = scp.ul.to_tree(df_out)
        df_out2 = scp.ul.to_cfmatrix(tree)
        df_out = df_out.loc[df_out2.index, df_out2.columns].copy()
        pd.testing.assert_frame_equal(df_out, df_out2, check_dtype=False)

    def test_tree_to_mtree(self):
        df_in = scp.datasets.test()
        df_out = scp.tl.phiscsb(df_in, alpha=0.0000001, beta=0.1)
        tree = scp.ul.to_tree(df_out)
        tree = scp.ul.to_mtree(tree)
        assert len(tree.nodes) == 13

    def test_to_tree(self, test_cf_data_1):
        data = scp.io.read(test_cf_data_1)
        tree = scp.ul.to_tree(data)
        assert len(list(tree.nodes)) == 10
        assert len(list(tree.edges)) == 9

    def test_mtree(self, test_cf_data_1):
        data = scp.io.read(test_cf_data_1)
        tree = scp.ul.to_tree(data)
        mtree = scp.ul.to_mtree(tree)
        assert len(mtree.nodes[8]["label"]) == 13

    @skip_graphviz
    def test_cells_muts_rooted_at(self, test_cf_data_1):
        data = scp.io.read(test_cf_data_1)

        tree = scp.ul.to_tree(data)
        assert len(tree.nodes) == 10
        scp.pl.clonal_tree(
            tree, show_id=True, muts_as_number=False, cells_as_number=False
        )

        tree = scp.pp.collapse(tree)
        assert len(tree.nodes) == 6
        scp.pl.clonal_tree(
            tree, show_id=True, muts_as_number=False, cells_as_number=False
        )

        res = scp.ul.cells_rooted_at(tree, "[8]")
        assert res.shape[0] == 33

        res = scp.ul.muts_rooted_at(tree, "[8]")
        assert res.shape[0] == 51

    def test_general(self, test_cf_data_1):
        with pytest.raises(RuntimeError):
            scp.ul.executable("kDPFC", "SPhyR")

        scp.ul.dir_base(test_cf_data_1)
        scp.ul.dirbase(test_cf_data_1)
        scp.ul.parse_params_file(
            "simNo_2-s_7-m_20-h_1-minVAF_0.1-ISAV_0-n_10-fp_0-fn_0.1-na_0-d_0-l_1000000"
            ".SC"
        )

        @scp.ul.timeit
        def _test():
            return None

        _test()

        assert True
