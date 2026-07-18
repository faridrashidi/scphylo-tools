"""Verify general-purpose matrix and tree utilities."""

import subprocess

import numpy as np
import pandas as pd
import pytest

import scphylo as scp


class TestUtils:
    """Exercise clustering, conversion, traversal, and helper utilities."""

    def test_dist_l1_ignore_na_read_only(self):
        """Verify L1 distance without mutating a read-only input matrix."""
        matrix = np.array([[0, 3, 1, 1], [1, 0, 1, 3]])
        expected = matrix.copy()
        matrix.setflags(write=False)

        dist = scp.ul.dist_l1_ignore_na(matrix)

        np.testing.assert_array_equal(dist, np.array([[0.0, 0.5], [0.5, 0.0]]))
        np.testing.assert_array_equal(matrix, expected)

    def test_hclustering_1(self):
        """Verify hierarchical clustering on a synthetic binary matrix."""
        np.random.seed(0)
        df = pd.DataFrame(np.random.randint(0, 2, size=(10, 10)))
        clusters = scp.ul.hclustering(df)
        assert clusters[6].value_counts()[5] == 3

    def test_hclustering_2(self):
        """Verify hierarchical clustering with the L1 metric."""
        adata = scp.datasets.example()
        clusters = scp.ul.hclustering(adata.to_df(), metric="l1")
        assert len(clusters) == 81

    def test_hclustering_3(self):
        """Verify hierarchical clustering with the cosine metric."""
        adata = scp.datasets.example()
        clusters = scp.ul.hclustering(adata.to_df(), metric="cosine")
        assert len(clusters) == 81

    def test_dist_dendro(self):
        """Verify the aggregate dendrogram distance for example data."""
        scp.settings.verbosity = 0
        adata = scp.datasets.example()
        dist = scp.ul.dist_dendro(adata)
        assert dist.sum() > 262456
        assert dist.sum() < 262457

    def test_tree_to_cfmatrix(self):
        """Verify a lossless tree-to-matrix round trip."""
        df_in = scp.datasets.test()
        df_out = scp.tl.phiscsb(df_in, alpha=0.0000001, beta=0.1)
        tree = scp.ul.to_tree(df_out)
        df_out2 = scp.ul.to_cfmatrix(tree)
        df_out = df_out.loc[df_out2.index, df_out2.columns].copy()
        pd.testing.assert_frame_equal(df_out, df_out2, check_dtype=False)

    def test_tree_to_mtree(self):
        """Verify conversion from a cell tree to a mutation tree."""
        df_in = scp.datasets.test()
        df_out = scp.tl.phiscsb(df_in, alpha=0.0000001, beta=0.1)
        tree = scp.ul.to_tree(df_out)
        tree = scp.ul.to_mtree(tree)
        assert len(tree.nodes) == 13

    def test_to_tree(self, test_cf_data_1):
        """Verify the node and edge counts of a converted tree."""
        data = scp.io.read(test_cf_data_1)
        tree = scp.ul.to_tree(data)
        assert len(list(tree.nodes)) == 10
        assert len(list(tree.edges)) == 9

    def test_mtree(self, test_cf_data_1):
        """Verify mutation labels after mutation-tree conversion."""
        data = scp.io.read(test_cf_data_1)
        tree = scp.ul.to_tree(data)
        mtree = scp.ul.to_mtree(tree)
        assert len(mtree.nodes[8]["label"]) == 13

    def test_cells_muts_rooted_at(self, test_cf_data_1):
        """Verify rooted cell and mutation queries after tree collapse."""
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
        """Verify executable lookup, path parsing, and timing helpers."""
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

    def test_external_tool_resolution_and_execution(self, tmp_path, monkeypatch):
        """Resolve and safely execute a tool whose path contains spaces."""
        tools_dir = tmp_path / "external tools"
        tools_dir.mkdir()
        executable = tools_dir / "fake-tool"
        executable.write_text(
            "#!/usr/bin/env python3\n"
            "import sys\n"
            "if sys.argv[1:] == ['fail']:\n"
            "    print('useful failure', file=sys.stderr)\n"
            "    raise SystemExit(7)\n"
            "print('|'.join(sys.argv[1:]))\n"
        )
        executable.chmod(0o755)

        resolved = scp.ul.resolve_executable(
            "fake-tool", "Fake Tool", tools_dir=tools_dir
        )
        result = scp.ul.run_external(
            [resolved, "argument with spaces"],
            "Fake Tool",
            stdout=subprocess.PIPE,
        )
        assert result.stdout.strip() == "argument with spaces"

        with pytest.raises(scp.ul.ExternalToolExecutionError, match="exit code 7"):
            scp.ul.run_external([resolved, "fail"], "Fake Tool")

        artifact = tools_dir / "tool.jar"
        artifact.write_text("not executable")
        artifact.chmod(0o600)
        assert scp.ul.resolve_external_file(
            "tool.jar", "Fake Tool", tools_dir=tools_dir
        ) == str(artifact.resolve())

        monkeypatch.setattr(scp.settings, "tools_dir", str(tools_dir))
        assert scp.ul.executable("fake-tool", "Fake Tool") == resolved
