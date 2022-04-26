"""
Visualizing a tree in clonal format
-----------------------------------

This example shows how to visualize an inferred tree.
"""

import scphylo as scp

# sphinx_gallery_thumbnail_path = "_static/thumbnails/clonal.png"

# %%
# First, we load a binary inferred single-cell genotype data.
inferred = scp.io.read(
    scp.ul.get_file("scphylo.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")
)

# %%
# Next we convert the inferred genotype matrix to a tree object.
tree = scp.ul.to_tree(inferred)

# %%
# Then we can draw the tree in `clonal` format i.e. mutations at the edges and cells at
# at the nodes of the tree.

# scp.pl.clonal_tree(tree)
# TODO: fix
