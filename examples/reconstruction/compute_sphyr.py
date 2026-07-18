"""SPhyR reconstruction.""" "\n---------------------"

# %%
# This example shows how to construct a phylogenetic tree using SPhyR on a binary
# single-cell genotype matrix.

import scphylo as scp

# sphinx_gallery_thumbnail_path = "_static/thumbnails/sphyr.png"

# %%
# First, we load a binary test single-cell genotype data.
df_in = scp.datasets.test()
df_in.head()

# %%
# Next, using :func:`scphylo.tl.sphyr` we remove the single-cell noises from the
# input. SPhyR requires CPLEX and ``sphyr_kDPFC``, supplied with ``executable``,
# ``SCPHYLO_TOOLS_DIR``, or ``PATH``. Regular documentation builds render this
# optional example without executing it; set ``SCPHYLO_RUN_EXTERNAL_EXAMPLES=sphyr``
# in a fully provisioned environment to run it (or use ``1`` to run every external
# example).
df_out = scp.tl.sphyr(df_in, alpha=0.0000001, beta=0.1)
df_out.head()


# %%
# Finally, using :func:`scphylo.ul.is_conflict_free_gusfield` we check whether the
# inferred genotype matrix is conflict-free or not.
is_cf = scp.ul.is_conflict_free_gusfield(df_out)
print(is_cf)
