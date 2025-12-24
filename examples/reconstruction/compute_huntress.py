"""
Reconstruct tree by HUNTRESS
----------------------------

This example shows how to construct a phylogenetic tree using HUNTRESS on a
binary single-cell genotype matrix.
"""

import scphylo as scp

# sphinx_gallery_thumbnail_path = "_static/thumbnails/huntress.png"

# %%
# First, we load a binary test single-cell genotype data.
df_in = scp.datasets.test()
df_in.head()

# %%
# Next, using :func:`scphylo.tl.huntress` we remove the single-cell noises from the
# input.
df_out = scp.tl.huntress(df_in, alpha=0.0000001, beta=0.1)
df_out.head()

# %%
# Finally, using :func:`scphylo.ul.is_conflict_free_gusfield` we check whether the
# inferred genotype matrix is conflict-free or not.
is_cf = scp.ul.is_conflict_free_gusfield(df_out)
print(is_cf)
