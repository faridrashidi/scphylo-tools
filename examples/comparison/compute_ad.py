"""Ancestor-descendent comparison.""" "\n-------------------------------"

# %%
# This example shows how to compare or measure two inferred genotype trees.

import scphylo as scp

# sphinx_gallery_thumbnail_path = "_static/thumbnails/ad.png"

# %%
# First, we load two binary test single-cell genotype data.
grnd = scp.io.read(
    scp.ul.get_file("scphylo.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")
)
sol = scp.io.read(
    scp.ul.get_file("scphylo.datasets/test/fp_1-fn_0.1-na_0.bnb.CFMatrix")
)

# %%
# Calculating the ancestor-descendent accuracy.
scp.tl.ad(grnd, sol)
