---
title: "scphylo-tools: A Python toolkit for single-cell tumor phylogenetic analysis"
tags:
  - python
  - biology
  - bioinformatics
  - cancer evolution
  - phylogenetics
  - single-cell sequencing
  - computational oncology
authors:
  - name: Farid Rashidi Mehrabadi
    orcid: 0000-0003-4103-4904
    affiliation: 1
affiliations:
  - name: Laboratory of Human Carcinogenesis, Center for Cancer Research, National Cancer Institute, National Institutes of Health, Bethesda, MD, USA
    index: 1
date: 11 January 2026
bibliography: paper.bib
---

# Summary

`scphylo-tools` is a Python library designed to unify single-cell tumor phylogeny inference methods. It addresses the lack of standardization in the field by providing a cohesive interface for data processing, tree reconstruction, visualization, and benchmarking. By streamlining these tasks, `scphylo-tools` empowers researchers with limited programming expertise to easily install and utilize complex inference methods, accelerates algorithm development, and ensures reproducible benchmarking of tumor phylogeny reconstruction methods.

# Statement of Need

Cancer is a dynamic evolutionary process driven by the acquisition of somatic mutations and the competitive selection of clonal populations. Single-Cell Sequencing (SCS) technologies have enabled the profiling of genomic alterations at cellular resolution, offering unprecedented insights into intratumor heterogeneity [@Navin_2011; @Wang_2014]. Elucidating the phylogenetic relationships between tumor cells is essential for understanding metastasis [@Leung_2017; @Roper_2020], drug resistance [@Kim_2014; @Gopalan_2021; @Gruen_2023], and clonal dynamics [@Laks_2019; @Liu_2025; @Hirsch_2025].

Reconstructing the evolutionary history of a tumor from SCS data presents unique computational challenges. These datasets are plagued by high noise levels, including Allele Drop-Out (ADO), false positives, missing data, and doublet artifacts [@ReviewBinary; @RashidiMehrabadi2022]. Consequently, a diverse array of computational tools has been developed to address these challenges, including stochastic methods like SCITE [@SCITE], infSCITE [@infSCITE], OncoNEM [@OncoNEM], SiFit [@SiFit], and SiCloneFit [@SiCloneFit]; combinatorial approaches such as PhISCS [@PhISCS], PhISCS-BnB [@PhISCS-BnB], SPhyR [@SPhyR], ScisTree [@ScisTree], gpps [@gpps], and SASC [@SASC]; and specialized methods including HUNTRESS [@HUNTRESS], SCIPhI [@SCIPhI], and Scelestial [@Scelestial].

However, the software landscape remains highly fragmented. Existing methods typically function as standalone binaries or scripts with inconsistent input/output formats, rendering comparative analysis difficult. Researchers attempting to utilize these tools face a laborious process of installation, data conversion, and script wrapping. Furthermore, integrating these tools into modern Python-based environments (e.g., alongside SCANPY [@SCANPY] or Biopython [@Biopython]) often requires custom development. `scphylo-tools` addresses these challenges by wrapping diverse state-of-the-art algorithms into a single, cohesive Python API.

The target audience includes computational biologists, bioinformaticians, and cancer researchers who need to:

- Compare multiple phylogeny inference methods on the same dataset
- Benchmark new algorithms against established methods
- Integrate tumor phylogeny analysis into existing Python workflows
- Access curated single-cell cancer datasets for research

# State of the Field

Several tools exist for tumor phylogeny inference, but none provide the unified framework that `scphylo-tools` offers. Individual tools like SCITE, PhISCS, and HUNTRESS each require separate installation, different input formats, and produce outputs that are not directly comparable. Workflow frameworks like Snakemake or Nextflow could orchestrate these tools, but would require significant custom development for format conversion and result standardization.

`scphylo-tools` distinguishes itself by:

1. **Unified Interface**: A consistent Python API across all wrapped methods, eliminating format conversion overhead
2. **Comprehensive Metrics**: Built-in implementation of specialized tumor tree comparison metrics (MLTD [@MLTD], CASet/DISC [@CASet_DISC], MP3 [@MP3])
3. **Curated Datasets**: Direct access to published SCS datasets from landmark cancer evolution studies [@Gawad_2014; @Morita_2020; @Hou_2012; @Leung_2017; @Wang_2014; @Wolf_2019]
4. **Ecosystem Integration**: Seamless integration with the Python scientific stack (NumPy, NetworkX, Matplotlib)

# Software Design

`scphylo-tools` follows a modular architecture inspired by established bioinformatics packages like SCANPY. The design philosophy prioritizes:

**Modular Organization**: The package is organized into functional modules (`io`, `pp`, `tl`, `pl`, `ul`, `datasets`) that mirror the phylogenetic analysis workflow. This separation allows users to engage only with the components relevant to their needs.

**Consistent Abstractions**: All solver wrappers expose a uniform interface, handling input formatting, binary execution, and output parsing internally. Trees are represented as NetworkX `DiGraph` objects, enabling interoperability with the broader Python ecosystem.

**Extensibility**: New inference methods can be added by implementing a standard wrapper interface, making the package straightforward to extend as new algorithms emerge.

**Reproducibility**: The `datasets` module provides versioned access to published datasets, and the simulation engine enables controlled generation of synthetic data with known ground truths for algorithm validation.

Key implementation decisions include using pandas DataFrames for genotype matrices (facilitating data manipulation), NetworkX for tree structures (enabling standard graph algorithms), and Matplotlib/Graphviz for visualization (providing publication-quality outputs).

# Research Impact Statement

`scphylo-tools` has been used in several published research studies:

- Analysis of melanoma subclonal evolution and immunotherapy resistance mechanisms [@Gruen_2023]
- Development of the Trisicell-PartF algorithm for evaluating inferred subclonal structures [@Trisicell-PartF]
- Investigation of expressed mutation profiles in single cells [@Trisicell-Boost]
- Metastatic migration pattern analysis using single-cell methylation sequencing [@Liu_2025; @Sgootr]
- Stochastic modeling of gene expression adaptation in tumor evolution [@Hirsch_2025]

The package is available on PyPI. Documentation and tutorials are hosted on Read the Docs, providing comprehensive guidance for new users.

# AI Usage Disclosure

No generative AI tools were used in the development of the `scphylo-tools` software, the writing of this manuscript, or the preparation of supporting documentation.

# Acknowledgements

We acknowledge the contributions of the open-source community and the authors of the underlying algorithms whose work is wrapped within this toolkit. This work was supported by the Intramural Research Program of the National Institutes of Health, National Cancer Institute.

# References
