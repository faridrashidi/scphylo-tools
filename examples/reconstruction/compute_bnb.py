"""
Reconstruct tree by PhISCS-BnB
------------------------------

This example shows how to construct a phylogenetic tree using PhISCS-BnB on a
binary single-cell genotype matrix.
"""

import scphylo as scp

# sphinx_gallery_thumbnail_path = "_static/thumbnails/bnb.png"

# %%
# First, we load a binary test single-cell genotype data.
df_in = scp.datasets.test()
df_in.head()

# %%
# Next, using :func:`scphylo.tl.bnb` we remove the single-cell noises from the
# input.

# TODO: fix
# df_out = scp.tl.bnb(df_in, bounding="simulated")
# df_out.head()


# %%
# Finally, using :func:`scphylo.ul.is_conflict_free_gusfield` we check whether the
# inferred genotype matrix is conflict-free or not.

# is_cf = scp.ul.is_conflict_free_gusfield(df_out)
# is_cf
