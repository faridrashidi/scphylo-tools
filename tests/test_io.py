import scphylo as scp


class TestIO:
    def test_read_newick(self, test_newick_1, test_newick_2):
        df1 = scp.io.read(test_newick_1)
        df2 = scp.io.read(test_newick_2)
        assert df1.shape == (185, 368)
        assert df2.shape == (185, 368)
