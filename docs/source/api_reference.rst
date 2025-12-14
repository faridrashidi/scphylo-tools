API
===

.. module:: scphylo

Import scphylo as::

   import scphylo as scp

After mutation calling and building the input data via our suggested
mutation calling pipeline.


Read/Write (io)
^^^^^^^^^^^^^^^
This module offers a bunch of functions for reading and writing of the data.

.. module:: scphylo.io
.. currentmodule:: scphylo
.. autosummary::
    :toctree: .

    io.read
    io.write


Preprocessing (pp)
^^^^^^^^^^^^^^^^^^
This module offers a bunch of functions for filtering and preprocessing of the
data.

.. module:: scphylo.pp
.. currentmodule:: scphylo
.. autosummary::
    :toctree: .

    pp.remove_mut_by_list
    pp.remove_cell_by_list
    pp.filter_mut_reference_must_present_in_at_least
    pp.filter_mut_mutant_must_present_in_at_least
    pp.bifiltering
    pp.consensus_combine


Tools (tl)
^^^^^^^^^^
This module offers a high-level API to compute the conflict-free solution
and calculating the probability of mutations seeding particular cells.

.. module:: scphylo.tl
.. currentmodule:: scphylo

**Solving the noisy input genotype matrix (scphylo-Boost)**

.. autosummary::
    :toctree: .

    tl.booster
    tl.scite
    tl.phiscsb
    tl.scistree
    tl.bnb
    tl.onconem
    tl.huntress
    tl.siclonefit
    tl.sphyr
    tl.gpps
    tl.phiscsi_bulk
    tl.scelestial


**Partition function calculation (scphylo-PartF)**

.. autosummary::
    :toctree: .

    tl.partition_function

**Consensus tree building (scphylo-Cons)**

.. autosummary::
    :toctree: .

    tl.consensus

**For comparing two phylogenetic trees**

.. autosummary::
    :toctree: .

    tl.ad
    tl.dl
    tl.mltd
    tl.tpted
    tl.caset
    tl.disc
    tl.mp3
    tl.rf
    tl.gs


Plotting (pl)
^^^^^^^^^^^^^
This module offers plotting the tree in clonal or dendrogram format.

.. module:: scphylo.pl
.. currentmodule:: scphylo
.. autosummary::
    :toctree: .

    pl.clonal_tree
    pl.dendro_tree


Utils (ul)
^^^^^^^^^^
This module offers a bunch of utility functions.

.. module:: scphylo.ul
.. currentmodule:: scphylo
.. autosummary::
    :toctree: .

    ul.to_tree
    ul.to_cfmatrix
    ul.to_mtree
    ul.hclustering
    ul.is_conflict_free
    ul.is_conflict_free_gusfield


Datasets (datasets)
^^^^^^^^^^^^^^^^^^^
This module offers a bunch of functions for simulating data.

.. module:: scphylo.datasets
.. currentmodule:: scphylo
.. autosummary::
    :toctree: .

    datasets.example
    datasets.simulate
    datasets.add_noise
    datasets.add_readcount
    datasets.melanoma20
    datasets.colorectal1
    datasets.colorectal2
    datasets.colorectal3
    datasets.acute_myeloid_leukemia1
    datasets.acute_myeloid_leukemia2
    datasets.acute_lymphocytic_leukemia1
    datasets.acute_lymphocytic_leukemia2
    datasets.acute_lymphocytic_leukemia3
    datasets.acute_lymphocytic_leukemia4
    datasets.acute_lymphocytic_leukemia5
    datasets.acute_lymphocytic_leukemia6
    datasets.high_grade_serous_ovarian_cancer1
    datasets.high_grade_serous_ovarian_cancer2
    datasets.high_grade_serous_ovarian_cancer3
    datasets.high_grade_serous_ovarian_cancer_3celllines
    datasets.myeloproliferative_neoplasms18
    datasets.myeloproliferative_neoplasms78
    datasets.myeloproliferative_neoplasms712
    datasets.renal_cell_carcinoma
    datasets.muscle_invasive_bladder
    datasets.erbc
    datasets.tnbc
