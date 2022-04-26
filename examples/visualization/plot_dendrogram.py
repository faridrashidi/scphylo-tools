"""
Visualizing a tree in dendrogram format
---------------------------------------

This example shows how to visualize an inferred tree.
"""

import scphylo as scp

# sphinx_gallery_thumbnail_path = "_static/thumbnails/dendrogram.png"

# %%
# First, we load a readcount single-cell genotype data and filter it.
adata = scp.datasets.example()
adata = adata[adata.obs.group.isin(["C16", "C11", "C22"]), :].copy()
scp.pp.filter_mut_vaf_greater_than_coverage_mutant_greater_than(
    adata, min_vaf=0.4, min_coverage_mutant=20, min_cells=2
)
scp.pp.filter_mut_reference_must_present_in_at_least(adata, min_cells=1)
scp.pp.filter_mut_mutant_must_present_in_at_least(adata, min_cells=2)
scp.pp.build_scmatrix(adata)
df_in = adata.to_df()

# %%
# Then we infer the tree using scistree algorithm.
df_out = scp.tl.scistree(df_in, alpha=0.001, beta=0.2)

# %%
# Next we convert the inferred genotype matrix to a tree object.
tree = scp.ul.to_tree(df_out)

# %%
# Finally we can draw the tree in `dendrogram` format.

# scp.pl.dendro_tree(
#     tree,
#     cell_info=adata.obs,
#     label_color="subclone_color",
#     width=1200,
#     height=600,
#     dpi=200,
#     distance_labels_to_bottom=3,
#     inner_node_type="both",
#     inner_node_size=2,
#     annotation=[
#         ("bar", "Axl", "Erbb3", 0.2),
#         ("bar", "Mitf", "Mitf", 0.2),
#     ],
# )
# TODO: fix
