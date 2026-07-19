import scphylo as scp


def example(is_expression=False):
    """Return an example for sanity checking and playing with scPhylo.

    is_expression : :obj:`bool`, optional
        Returns the expression dataset instead of the genotype one, by default False

    Returns
    -------
    :class:`anndata.AnnData`
        An object that cells are in `.obs` and mutations are in `.var`.
    """
    if is_expression:
        return scp.io.read(scp.ul.get_file("scphylo.datasets/data/expression.h5ad.gz"))
    else:
        return scp.io.read(scp.ul.get_file("scphylo.datasets/data/genotype.h5ad.gz"))


def test():
    """Load the small genotype dataset used by the package tests."""
    df = scp.io.read(scp.ul.get_file("scphylo.datasets/test/test.tsv"))
    return df


def melanoma20():
    """Mouse Melanoma dataset with 20 sublines.

    This dataset was introduced in :cite:`Wolf_2019` and was used in:

    * :cite:`PhISCS-BnB` Figure 1.

    The size is n_cells × n_muts = 20 × 2367

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.

        `.layers['solution_fig1']` is the solution presented in Figure 1 of PhISCS-BnB
        paper.
    """
    adata = scp.io.read(scp.ul.get_file("scphylo.datasets/real/melanoma20.h5ad"))
    return adata


def colorectal1():
    """Human Colorectal Cancer (Patient 1).

    This dataset was introduced in :cite:`Leung_2017` and was used in:

    * :cite:`B-SCITE` Figure 8a.
    * :cite:`SiFit` Figure 6.
    * :cite:`SPhyR` Table 1.
    * :cite:`SiCloneFit` Figure 3.

    The size is n_cells × n_muts = 178 × 16

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.

    Notes
    -----
    This dataset includes single cells from two sites of the patient body;
    133 single cells from colon as primary tumor site and 45 single cells from liver
    as the tumor metastatic site (178 in total). The number of mutations in this
    dataset is 16. One can remove the cells in this dataset that carry none of these 16
    mutations before feeding it to our network. After removing cells with zero profile,
    the number of cells are 40 and 32 from primary and metastatic sites, respectively
    (72 in total).
    """
    adata = scp.io.read(scp.ul.get_file("scphylo.datasets/real/colorectal1.h5ad"))
    return adata


def colorectal2(readcount=False):
    """Human Colorectal Cancer (Patient 2).

    This dataset was introduced in :cite:`Leung_2017` and was used in:

    * :cite:`PhISCS` Figure 7.
    * :cite:`B-SCITE` Figure 8b.
    * :cite:`SiCloneFit` Figure 4.
    * :cite:`SCARLET` Figure 4.

    The dataset has size n_cells × n_muts = 182 × 36.

    Parameters
    ----------
    readcount : :obj:`bool`, optional
        Include mutant and total single-cell read-count layers.

    Returns
    -------
    :class:`anndata.AnnData`
        An AnnData object in which `.X` contains the observed genotype calls: 0 is
        wild type, 1 is mutant, and 3 is missing.

            - `.layers['mutant']` and `.layers['total']` contain read counts when
                ``readcount=True``.
            - `.obs['copy_number_profile']` contains the SCARLET copy-number profile.
            - `.obs['phiscs_fig7']` and `.var['phiscs_fig7']` identify the historical
                78 × 25 PhISCS subset.
            - `.uns['phiscs_fig7']` contains that subset's cell and mutation names,
                Figure 7 solutions, and parameters.
            - `.var` contains genomic coordinates and available bulk-sample counts.
            - `.uns['provenance']` records the genotype and read-count sources.

    Notes
    -----
    The full genotype matrix is reconstructed from Supplementary Figure 7 of
    :cite:`Leung_2017`. The historical 78 × 25 matrix and PhISCS solutions remain
    available as the name-indexed ``.uns['phiscs_fig7']`` derivative. Bulk counts
    are unavailable for the 11 loci excluded from that derivative and are stored as
    missing values in ``.var``.

    The read counts come from the processed SCARLET data. That source labels SRA run
    SRR3472637 as MA-84, whereas the Leung tables and Figure 7 label the same sample
    MA_94. The source-only low-coverage cells PA_31, PDD_21, PDD_29, and PDD_66 are
    excluded to match the 182 post-QC Figure 7 cells.
    """
    if readcount:
        adata = scp.io.read(
            scp.ul.get_file("scphylo.datasets/real/colorectal2.rc.h5ad")
        )
    else:
        adata = scp.io.read(scp.ul.get_file("scphylo.datasets/real/colorectal2.h5ad"))
    return adata


def colorectal3():
    """Human Colorectal Cancer (Patient CRC0827).

    This dataset was introduced in :cite:`Wu_2016` and was used in:

    * :cite:`SiFit` Figure 5.

    The complete published matrix has size n_cells × n_muts = 61 × 77.

    Returns
    -------
    :class:`anndata.AnnData`
        An AnnData object in which ``.X`` contains the observed genotype calls:
        0 is reference, 1 is mutation present, and 3 is missing/no coverage.

            - ``.obs['source_tissue']`` distinguishes the 35 cancer-biopsy,
                13 adenomatous-polyp, and 13 matched-normal colorectal cells.
            - ``.obs['non_normal_tissue']`` identifies the historical 48 × 77
                view obtained by excluding the 13 matched-normal cells.
            - ``.var['sifit_index']`` preserves the mutation numbering in
                Supplementary Figure S4a.
            - ``.uns['provenance']`` records the source figures, exact raster
                extraction, validation, and archival SiFit code.

    Notes
    -----
    The previous 48 × 77 documentation described the non-normal-tissue slice,
    not the input to :cite:`SiFit` Figure 5. Both the paper and its archived
    CRC0827 runner use all 61 cells: 35 cancer-biopsy cells, 13
    adenomatous-polyp cells, and 13 matched-normal controls. The complete matrix
    is therefore returned by default, while the 48-cell view remains available
    through the name-indexed ``non_normal_tissue`` mask.

    The calls were decoded from the lossless raster embedded in :cite:`SiFit`
    Supplementary Figure S4a. The extraction reproduces the paper's 9.4% missing
    rate and 1,847 four-gamete violations among all 2,926 mutation pairs.
    """
    adata = scp.io.read(scp.ul.get_file("scphylo.datasets/real/colorectal3.h5ad"))
    return adata


def acute_lymphocytic_leukemia1():
    """Human Acute Lymphocytic Leukemia dataset (Patient 1).

    This dataset was introduced in :cite:`Gawad_2014` and was used in:

    * :cite:`B-SCITE` Figure 5.
    * :cite:`infSCITE` Figure S16.

    The size is n_cells × n_muts = 111 × 20

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.
    """
    adata = scp.io.read(
        scp.ul.get_file("scphylo.datasets/real/acute_lymphocytic_leukemia1.h5ad")
    )
    return adata


def acute_lymphocytic_leukemia2():
    """Human Acute Lymphocytic Leukemia dataset (Patient 2).

    This dataset was introduced in :cite:`Gawad_2014` and was used in:

    * :cite:`PhISCS` Figure 9.
    * :cite:`B-SCITE` in Figure 6.
    * :cite:`infSCITE` Figure S17.
    * :cite:`Phyolin` Table 2.

    The dataset has size n_cells × n_muts = 115 × 16.

    Returns
    -------
    :class:`anndata.AnnData`
        An AnnData object in which `.X` contains the complete noisy genotype matrix.

            - `.obs['phiscs_fig9']` and `.var['phiscs_fig9']` identify the historical
                102 × 16 PhISCS subset.
            - `.uns['phiscs_fig9']` contains that subset's cell and mutation names,
                Figure 9 solution, and parameters.
            - `.var` contains genomic coordinates and bulk-sample counts.
            - `.uns['provenance']` records the genotype and PhISCS preprocessing
                sources.

    Notes
    -----
    The complete 115-cell input was analyzed in :cite:`B-SCITE` Figure 6 and
    :cite:`infSCITE` Figure S17. Before running PhISCS, Single Cell Genotyper
    identified and removed 13 predicted doublets, leaving the 102 cells used in
    :cite:`PhISCS` Figure 9. That name-indexed input, solution, and parameter set
    remain available in ``.uns['phiscs_fig9']`` without padding the published result.
    """
    adata = scp.io.read(
        scp.ul.get_file("scphylo.datasets/real/acute_lymphocytic_leukemia2.h5ad")
    )
    return adata


def acute_lymphocytic_leukemia3():
    """Human Acute Lymphocytic Leukemia dataset (Patient 3).

    The bundled binary genotype matrix was introduced in :cite:`Gawad_2014`
    and was used in:

    * :cite:`infSCITE` Figure S18.

    The dataset has size n_cells × n_muts = 150 × 49.

    Returns
    -------
    :class:`anndata.AnnData`
        An AnnData object in which `.X` contains the complete published noisy
        genotype matrix.

            - `.obs['gawad_cluster']` contains the source's five cluster labels.
            - `.var` contains genomic coordinates and bulk-sample counts.
            - `.uns['provenance']` records the genotype source and the later
                255-cell raw-read reanalysis.

    Notes
    -----
    :cite:`Gawad_2014` retained the 150 Patient 3 cells whose estimated allelic
    dropout was below 30%. :cite:`SCIPhI` Figure 5 and :cite:`ScisTree` Figure S3
    instead reanalyzed 255 Patient 3 raw-read libraries from accession SRP044380
    and independently called variants. Those data are not 105 omitted rows of
    this 49-mutation matrix; ScisTree, for example, reports 406 called SNV sites.
    A 255-cell read-count or probability dataset would therefore be a separate
    dataset, not a dimension correction to this one.
    """
    adata = scp.io.read(
        scp.ul.get_file("scphylo.datasets/real/acute_lymphocytic_leukemia3.h5ad")
    )
    return adata


def acute_lymphocytic_leukemia4():
    """Human Acute Lymphocytic Leukemia dataset (Patient 4).

    This dataset was introduced in :cite:`Gawad_2014` and was used in:

    * :cite:`infSCITE` Figure S19.
    * :cite:`gpps` Figure 3.
    * :cite:`SASC` Figure 6.

    The size is n_cells × n_muts = 143 × 78

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.
    """
    adata = scp.io.read(
        scp.ul.get_file("scphylo.datasets/real/acute_lymphocytic_leukemia4.h5ad")
    )
    return adata


def acute_lymphocytic_leukemia5():
    """Human Acute Lymphocytic Leukemia dataset (Patient 5).

    This dataset was introduced in :cite:`Gawad_2014` and was used in:

    * :cite:`infSCITE` Figure S20.
    * :cite:`SASC` Figure 7.

    The size is n_cells × n_muts = 96 × 105

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.
    """
    adata = scp.io.read(
        scp.ul.get_file("scphylo.datasets/real/acute_lymphocytic_leukemia5.h5ad")
    )
    return adata


def acute_lymphocytic_leukemia6():
    """Human Acute Lymphocytic Leukemia dataset (Patient 6).

    This dataset was introduced in :cite:`Gawad_2014` and was used in:

    * :cite:`infSCITE` Figure S21.
    * :cite:`Phyolin` Table 2.

    The size is n_cells × n_muts = 146 × 10

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.
    """
    adata = scp.io.read(
        scp.ul.get_file("scphylo.datasets/real/acute_lymphocytic_leukemia6.h5ad")
    )
    return adata


def tnbc():
    """Triple-negative Breast Cancer.

    This dataset was introduced in :cite:`Wang_2014` and was used in:

    * :cite:`SCIPhI` Figure 4.
    * :cite:`B-SCITE` Figure 7.
    * :cite:`TRaIT` Figure 6.

    The size is n_cells × n_muts = 16 × 20

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.

    Examples
    --------
    >>> adata = scp.datasets.tnbc()
    >>> df_in = adata.to_df()
    """
    adata = scp.io.read(scp.ul.get_file("scphylo.datasets/real/tnbc.h5ad"))
    return adata


def erbc():
    """Oestrogen-receptor-positive (ER+) Breast Cancer.

    This dataset was introduced in :cite:`Wang_2014` and was used in:

    * :cite:`SCITE` Figure S8 and S9.
    * :cite:`infSCITE` Figure S15.
    * :cite:`gpps` Figure 1.

    The size is n_cells × n_muts = 47 × 40

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.

            - `.uns['params_scite']` is parameters inferred by SCITE.

    Examples
    --------
    >>> adata = scp.datasets.erbc()
    >>> df_in = adata.to_df()
    """
    adata = scp.io.read(scp.ul.get_file("scphylo.datasets/real/erbc.h5ad"))
    return adata


def muscle_invasive_bladder():
    """Muscle Invasive Bladder Cancer.

    This dataset was introduced in :cite:`Li_2012` and was used in:

    * :cite:`OncoNEM` Figure 6B.

    The size is n_cells × n_muts = 44 × 443

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.

            - `.uns['params_onconem']` is parameters inferred by OncoNEM.

    Examples
    --------
    >>> adata = scp.datasets.muscle_invasive_bladder()
    >>> df_in = adata.to_df()
    """
    adata = scp.io.read(
        scp.ul.get_file("scphylo.datasets/real/muscle_invasive_bladder.h5ad")
    )
    return adata


def renal_cell_carcinoma():
    """Clear-cell Renal-cell Carcinoma.

    This dataset was introduced in :cite:`Xu_2012` and was used in:

    * :cite:`SCITE` Figure S6 and S7.
    * :cite:`infSCITE` Figure S14.

    The size is n_cells × n_muts = 17 × 35

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.

            - `.uns['params_scite']` is parameters inferred by SCITE.

    Examples
    --------
    >>> adata = scp.datasets.renal_cell_carcinoma()
    >>> df_in = adata.to_df()
    """
    adata = scp.io.read(
        scp.ul.get_file("scphylo.datasets/real/renal_cell_carcinoma.h5ad")
    )
    return adata


def myeloproliferative_neoplasms712():
    """JAK2-Negative Myeloproliferative Neoplasm.

    This dataset was introduced in :cite:`Hou_2012` and was used in:

    * :cite:`OncoNEM` Figure 6D.

    The size is n_cells × n_muts = 58 × 712

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.

            - `.uns['params_onconem']` is parameters inferred by OncoNEM.
    """
    adata = scp.io.read(
        scp.ul.get_file("scphylo.datasets/real/myeloproliferative_neoplasms712.h5ad")
    )
    return adata


def myeloproliferative_neoplasms78():
    """JAK2-Negative Myeloproliferative Neoplasm.

    This dataset was introduced in :cite:`Hou_2012` and was used in:

    * :cite:`SCITE` Figure S5.

    The size is n_cells × n_muts = 58 × 78

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.

    Notes
    -----
    The original dataset contains 712 mutations but 78 ones were considered as
    non-synonymous mutations from the full data.
    """
    adata = scp.io.read(
        scp.ul.get_file("scphylo.datasets/real/myeloproliferative_neoplasms78.h5ad")
    )
    return adata


def myeloproliferative_neoplasms18():
    """JAK2-Negative Myeloproliferative Neoplasm.

    This dataset was introduced in :cite:`Hou_2012` and was used in:

    * :cite:`SCITE` Figure S2, S3 and S4.
    * :cite:`Kim_2014` Figure 1.
    * :cite:`infSCITE` Figure S13.
    * :cite:`gpps` Figure 2.

    The size is n_cells × n_muts = 58 × 18

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.

            - `.uns['params_scite']` is parameters inferred by SCITE.

    Notes
    -----
    The original dataset contains 712 mutations but 18 ones were considered as
    cancer related mutations from the full data.
    """
    adata = scp.io.read(
        scp.ul.get_file("scphylo.datasets/real/myeloproliferative_neoplasms18.h5ad")
    )
    return adata


def high_grade_serous_ovarian_cancer1():
    """High Grade Serous Ovarian Cancer (Patient 2).

    This dataset was introduced in :cite:`McPherson_2016` and was used in:

    * :cite:`infSCITE` Figure S22.

    The size is n_cells × n_muts = 588 × 37

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.
    """
    # adata = scp.io.read(
    #     scp.ul.get_file(
    #         "scphylo.datasets/real/high_grade_serous_ovarian_cancer1.h5ad"
    #     )
    # )
    # TODO: extract
    return None


def high_grade_serous_ovarian_cancer2():
    """High Grade Serous Ovarian Cancer (Patient 3).

    This dataset was introduced in :cite:`McPherson_2016` and was used in:

    * :cite:`infSCITE` Figure S23.

    The size is n_cells × n_muts = 672 × 60

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.
    """
    # adata = scp.io.read(
    #     scp.ul.get_file(
    #         "scphylo.datasets/real/high_grade_serous_ovarian_cancer2.h5ad"
    #     )
    # )
    # TODO: extract
    return None


def high_grade_serous_ovarian_cancer3():
    """High Grade Serous Ovarian Cancer (Patient 9).

    This dataset was introduced in :cite:`McPherson_2016` and was used in:

    * :cite:`infSCITE` Figure S24.
    * :cite:`scVILP` Figure 5.
    * :cite:`SCIPhI` Figure S10.

    The complete targeted single-nucleotide-locus assay has size
    n_cells × n_loci = 420 × 43.

    Returns
    -------
    :class:`anndata.AnnData`
        An AnnData object in which ``.X`` is a binary alternate-allele-presence
        projection of the published allele calls: 0 is reference-only (A), 1 means
        the alternate allele is present (AB or B), and 3 is uncalled or missing.

            - ``.layers['allele_state']`` preserves the original calls as 0=A,
                1=AB, 2=B, and 3=missing.
            - ``.layers['mutant']`` and ``.layers['total']`` contain
                alternate-supporting and total source read counts, respectively.
            - ``.layers['ref_p_value']`` and ``.layers['alt_p_value']`` contain the
                source binomial exact-test p-values for allele presence.
            - ``.obs['mcpherson_clonal_analysis']`` identifies the 402 nuclei used
                in the original clonal analysis.
            - ``.obs['sciphi_scvilp']`` identifies the later 370-cell read-depth
                subset.
            - ``.var['infscite_fig_s24']`` identifies the 37 somatic-SNV targets
                used in the primary infSCITE analysis.
            - ``.uns['provenance']`` records the source and derivative filters.

    Notes
    -----
    The apparent published dimensions describe different, reproducible views of
    the same assay. The source contains 420 nuclei and 43 targeted loci: 37 somatic
    SNVs and six germline heterozygous normal-marker SNPs putatively lost in tumor
    cells. :cite:`infSCITE` Figure S24 excluded the six normal markers, giving
    420 × 37. :cite:`SCIPhI` Figure S10 retained all 43 loci but selected nuclei
    with more than 10,000 reads summed across them, giving 370 × 43;
    :cite:`scVILP` Figure 5 uses the same read-count matrix. Both derivatives are
    retained as name-indexed masks rather than replacing the full source panel.
    """
    adata = scp.io.read(
        scp.ul.get_file("scphylo.datasets/real/high_grade_serous_ovarian_cancer3.h5ad")
    )
    return adata


def high_grade_serous_ovarian_cancer_3celllines():
    """High Grade Serous Ovarian Cancer (3 cell lines).

    This dataset was introduced in :cite:`Laks_2019`
    and preprocessed in :cite:`McPherson_2019`.
    The phylogeny is presented in Figure 3H.

    The size is n_cells × n_muts = 891 × 13666

    Note that 402 mutations were deleted during evolution
    in the inferred tree by original study So, here they were filtered out
    (meaning they are 14068 mutations in total in the original study).

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy
        (it is obtained based on the number of mutant/total reads).

            - `.obs['clone_id']` is the clone id to which the cell is assigned in
                Figure 3H.
            - `.obs['group_color']` is unique colors for each 'clone_id'.
            - `.obs['cell_name']` is a new name for each cell based on the
                'group_color'.
            - `.layers['mutant']` is the number of mutant reads at each locus in each
                cell.
            - `.layers['total']` is the total number of reads at each locus in each
                cell.
            - `.layers['ground']` is the solution inferred in Figure 3H of the original
                paper.
            - `.uns['params_ground']` is parameters inferred by comparing ground and
                noisy matrices.
            - `.var` includes information of the bulk samples.

    Examples
    --------
    >>> sc = scp.datasets.high_grade_serous_ovarian_cancer_3celllines()
    >>> print(sc)
    AnnData object with n_obs × n_vars = 891 × 13666
        obs: 'clone_id', 'group_color', 'cell_name'
        uns: 'params_ground'
        layers: 'ground', 'mutant', 'total'
    """
    adata = scp.io.read(scp.ul.get_file("scphylo.datasets/real/ovarian.h5ad.gz"))
    adata.var_names = adata.var_names.str.replace(":", "_")
    return adata


def oligodendroglioma_idh_mutated_tumor():
    """Oligodendroglioma IDH-mutated tumor.

    This dataset was introduced in :cite:`Tirosh_2016` and was used in:

    * :cite:`SASC` Figure 5.

    The size is n_cells × n_muts = 579 × 77

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.
    """
    adata = scp.io.read(
        scp.ul.get_file(
            "scphylo.datasets/real/oligodendroglioma_idh_mutated_tumor.h5ad"
        )
    )
    return adata


def acute_lymphocytic_leukemia_many():
    """Human Acute Lymphocytic Leukemia datasets.

    This dataset was introduced in :cite:`Morita_2020` and was used in:

    * :cite:`Phyolin` Figure 3 and Table 1.

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.
    """
    # TODO: extract all
    return None


def isogenic_fibroblast_cell_line():
    """Isogenic Fibroblast cell line dataset.

    This dataset was introduced in :cite:`Leung_2015` and was used in:

    * :cite:`SCIPhI` Figure S9.

    The size is n_cells × n_muts = 19 × ?

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.
    """
    # TODO: extract
    return None


def acute_myeloid_leukemia1():
    """Human Acute Myeloid Leukemia dataset (Patient 1).

    This dataset was introduced in :cite:`Velten_2021` and was used in:

    The size is n_cells × n_muts = 1430 × 15

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.
    """
    adata = scp.io.read(
        scp.ul.get_file("scphylo.datasets/real/acute_myeloid_leukemia1.h5ad")
    )
    return adata


def acute_myeloid_leukemia2():
    """Human Acute Myeloid Leukemia dataset (Patient 2).

    This dataset was introduced in :cite:`Velten_2021` and was used in:

    The size is n_cells × n_muts = 1066 × 21

    Returns
    -------
    :class:`anndata.AnnData`
        An anndata in which `.X` is the input noisy.
    """
    adata = scp.io.read(
        scp.ul.get_file("scphylo.datasets/real/acute_myeloid_leukemia2.h5ad")
    )
    return adata
