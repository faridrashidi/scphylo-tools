"""
Reconstruct tree by Trisicell-Boost
-----------------------------------

This example shows how to construct a phylogenetic tree using Trisicell-Boost on a
binary single-cell genotype matrix.
"""

import scphylo as scp

# sphinx_gallery_thumbnail_path = "_static/thumbnails/trisicell-boost.png"

# %%
# First, we load a binary test single-cell genotype data.
df_in = scp.datasets.test()
df_in.head()

# %%
# Next, using :func:`scphylo.tl.booster` we remove the single-cell noises from the
# input.
df_out = scp.tl.booster(
    df_in,
    alpha=0.0000001,
    beta=0.1,
    solver="PhISCS",
    sample_on="muts",
    sample_size=15,
    n_samples=88,
    begin_index=0,
    n_jobs=1,
    dep_weight=50,
)
df_out.head()

# %%
# Finally, using :func:`scphylo.ul.is_conflict_free_gusfield` we check whether the
# inferred genotype matrix is conflict-free or not.
is_cf = scp.ul.is_conflict_free_gusfield(df_out)
print(is_cf)
