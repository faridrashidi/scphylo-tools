"""
Comparing scores for two phylogenetic trees
-------------------------------------------

This example shows how to compare/measure two inferred genotype data (trees).
"""

import scphylo as scp

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

# %%
# Calculating the different-lineage accuracy.
scp.tl.dl(grnd, sol)

# %%
# Calculating the multi-labeled tree dissimilarity measure (MLTD).
scp.tl.mltd(grnd, sol)

# %%
# Calculating the tumor phylogeny tree edit distance measure (TPTED).
scp.tl.tpted(grnd, sol)

# %%
# Calculating the distinctly inherited sets score (DISC).
scp.tl.disc(grnd, sol)

# %%
# Calculating the commonly ancestor sets score (CASet).
scp.tl.caset(grnd, sol)

# %%
# Calculating the Triplet-based similarity score (MP3).
scp.tl.mp3(grnd, sol)

# %%
# Calculating the Robinsold-Foulds similarity score (1 - normalized_distance).
scp.tl.rf(grnd, sol)
