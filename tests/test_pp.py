import scphylo as scp

from ._helpers import skip_gurobi


class TestPreProcessing:
    @skip_gurobi
    def test_bifiltering(self):
        df_in = scp.datasets.test()
        df_filtered = scp.pp.bifiltering(df_in, 0.5, 0.2)
        assert df_filtered.shape == (10, 4)

    def test_binary(self):
        scp.settings.verbosity = (
            0  # UnicodeEncodeError: 'ascii' codec can't encode character '\xd7'
        )
        df_in = scp.datasets.test()
        scp.pp.binarym_filter_private_mutations(df_in)
        scp.pp.binarym_filter_clonal_mutations(df_in)
        scp.pp.binarym_filter_nonsense_mutations(df_in)
        scp.pp.binarym_statistics(df_in)
        assert df_in.shape == (20, 19)
        adata = scp.datasets.example()
        scp.pp.build_scmatrix(adata)
        df_in = scp.pp.consensus_combine(adata.to_df())
        assert df_in.shape == (11, 452)

    def test_readcount(self):
        scp.settings.verbosity = (
            0  # UnicodeEncodeError: 'ascii' codec can't encode character '\xd7'
        )
        adata = scp.datasets.example()
        scp.pp.filter_mut_vaf_greater_than_coverage_mutant_greater_than(
            adata, min_vaf=0.4, min_coverage_mutant=20, min_cells=2
        )
        scp.pp.filter_mut_reference_must_present_in_at_least(adata, min_cells=1)
        scp.pp.filter_mut_mutant_must_present_in_at_least(adata, min_cells=2)
        scp.pp.statistics(adata)
        assert adata.shape == (83, 267)
        scp.pp.group_obs_apply_func(adata, group_key="group")
        scp.pp.remove_cell_by_list(adata, ["C15_1"])
        scp.pp.keep_cell_by_list(adata, ["C15_2", "C15_3"])
        assert adata.shape == (2, 267)

    def test_tree(self):
        adata = scp.datasets.high_grade_serous_ovarian_cancer_3celllines()
        df_out = adata.to_df(layer="ground")[adata.var_names[:1000]]
        tree = scp.ul.to_tree(df_out)
        tree = scp.pp.collapse(tree)
        sampled_cells = scp.pp.sample_from_tree(tree, 0.2, axis="cell")
        assert len(sampled_cells) > 160
