"""
Reconstruct tree by Trisicell-PartF
-----------------------------------

This example shows how to assess a phylogenetic tree using Trisicell-PartF on a
binary single-cell genotype matrix.
"""

import scphylo as scp

# sphinx_gallery_thumbnail_path = "_static/thumbnails/trisicell-partf.png"

# %%
# First, we load a binary test single-cell genotype data.
df_in = scp.datasets.test()
df_in.head()

# %%
# Then, using :func:`scphylo.tl.partf` we calculate the probability of `cell6` and
# `cell17` seeded by `mut12`
probs = scp.tl.partition_function(
    df_in,
    alpha=0.000001,
    beta=0.1,
    n_samples=100,
    n_batches=10,
    muts=["mut12"],
    cells=["cell6", "cell17"],
)
probs.mean(axis=1).round(4).values[0]
