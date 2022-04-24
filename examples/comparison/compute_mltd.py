"""
Comparing two trees by MLTD
---------------------------

This example shows how to compare/measure two inferred genotype data (trees).
"""

import scphylo as scp

# sphinx_gallery_thumbnail_path = "_static/thumbnails/mltd.png"

# %%
# First, we load two binary test single-cell genotype data.
grnd = scp.io.read(
    scp.ul.get_file("scphylo.datasets/test/fp_0-fn_0-na_0.ground.CFMatrix")
)
sol = scp.io.read(
    scp.ul.get_file("scphylo.datasets/test/fp_1-fn_0.1-na_0.bnb.CFMatrix")
)

# %%
# Calculating the multi-labeled tree dissimilarity measure (MLTD).
scp.tl.mltd(grnd, sol)
