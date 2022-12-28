import scphylo as scp

from ._helpers import skip_graphviz, skip_rpy2


class TestPlotting:
    def test_heatmap(self):
        adata = scp.datasets.example()
        adata = adata[adata.obs.group.isin(["C16", "C11", "C22"]), :].copy()
        scp.pp.filter_mut_vaf_greater_than_coverage_mutant_greater_than(
            adata, min_vaf=0.4, min_coverage_mutant=20, min_cells=2
        )
        scp.pp.filter_mut_reference_must_present_in_at_least(adata, min_cells=1)
        scp.pp.filter_mut_mutant_must_present_in_at_least(adata, min_cells=2)
        scp.pp.build_scmatrix(adata)
        scp.pl.heatmap(adata, color_attrs="subclone_color")
        assert True
        df_in = adata.to_df()
        df_out = scp.tl.scistree(df_in, alpha=0.001, beta=0.2)
        D = 1 * (df_in == 1) * (df_out == 0) + 2 * (df_in == 0) * (df_out == 1)
        adata.obsm["flips"] = D
        scp.pl.heatmap(
            adata,
            color_attrs="subclone_color",
            layer="flips",
            rvb=["#FFFFFF", "#D92347", "#A9D0F5"],
            vmin=0,
            vmax=2,
        )
        assert True

    @skip_graphviz
    def test_clonal_tree(self, test_cf_data_1):
        data = scp.io.read(test_cf_data_1)
        tree = scp.ul.to_tree(data)
        scp.pl.clonal_tree(tree)

    @skip_graphviz
    def test_clonal_tree_with_coloring(self):
        adata = scp.datasets.high_grade_serous_ovarian_cancer_3celllines()
        df_out = adata.to_df(layer="ground")[adata.var_names[:1000]]
        tree = scp.ul.to_tree(df_out)
        scp.pl.clonal_tree(
            tree,
            muts_as_number=True,
            cells_as_number=False,
            cell_info=adata.obs,
            color_attr="group_color",
        )
        scp.pl.clonal_tree(
            tree,
            muts_as_number=True,
            cells_as_number=True,
            show_id=True,
        )
        assert True

    @skip_rpy2()
    def test_dendro_tree_1(self, test_cf_data_1):
        data = scp.io.read(test_cf_data_1)
        tree = scp.ul.to_tree(data)
        scp.pl.dendro_tree(tree)
        assert True

    @skip_rpy2()
    def test_dendro_tree_2(self):
        adata = scp.datasets.example()
        adata = adata[adata.obs.group.isin(["C16", "C11", "C22"]), :].copy()
        scp.pp.filter_mut_vaf_greater_than_coverage_mutant_greater_than(
            adata, min_vaf=0.4, min_coverage_mutant=20, min_cells=2
        )
        scp.pp.filter_mut_reference_must_present_in_at_least(adata, min_cells=1)
        scp.pp.filter_mut_mutant_must_present_in_at_least(adata, min_cells=2)
        scp.pp.build_scmatrix(adata)
        df_in = adata.to_df()
        df_out = scp.tl.scistree(df_in, alpha=0.001, beta=0.2)
        tree = scp.ul.to_tree(df_out)
        scp.pl.dendro_tree(
            tree,
            cell_info=adata.obs,
            label_color="subclone_color",
            width=1200,
            height=600,
            dpi=200,
            distance_labels_to_bottom=3,
            inner_node_type="both",
            inner_node_size=2,
            annotation=[
                ("bar", "Axl", "Erbb3", 0.2),
                ("bar", "Mitf", "Mitf", 0.2),
            ],
        )
        assert True
