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
  - name: Laboratory of Human Carcinogenesis, Center for Cancer Research, National Cancer Institute, National Institutes of Health, Bethesda, MD, United States of America
    index: 1
date: 11 January 2026
bibliography: paper.bib
---

# Summary

`scphylo-tools` is a Python library designed to unify single-cell tumor phylogeny inference methods. It addresses the lack of standardization in the field by providing a cohesive interface for data processing, tree reconstruction, visualization, and benchmarking. By streamlining these tasks, `scphylo-tools` empowers researchers with limited programming expertise to easily install and utilize complex inference methods, accelerates algorithm development, and ensures reproducible benchmarking of tumor phylogeny reconstruction methods.

# Statement of Need

Cancer is a dynamic evolutionary process driven by the acquisition of somatic mutations and the competitive selection of clonal populations. Single-Cell Sequencing (SCS) technologies have enabled the profiling of genomic alterations at cellular resolution, offering unprecedented insights into intratumor heterogeneity [@Navin_2011; @Wang_2014]. Elucidating the phylogenetic relationships between tumor cells is essential for understanding metastasis [@Leung_2017; @Roper_2020], drug resistance [@Kim_2014; @Gopalan_2021; @Gruen_2023], and clonal dynamics [@Laks_2019; @Liu_2025; @Hirsch_2025].

Reconstructing the evolutionary history of a tumor from SCS data presents unique computational challenges. These datasets are plagued by high noise levels, including Allele Drop-Out (ADO), false positives, missing data, and doublet artifacts [@ReviewBinary; @RashidiMehrabadi2022]. Consequently, a diverse array of computational tools has been developed to address these challenges, including stochastic methods like SCITE [@SCITE], infSCITE [@infSCITE], OncoNEM [@OncoNEM], SiFit [@SiFit], and SiCloneFit [@SiCloneFit]; combinatorial approaches such as PhISCS [@PhISCS], PhISCS-BnB [@PhISCS-BnB], SPhyR [@SPhyR], ScisTree [@ScisTree], gpps [@gpps], and SASC [@SASC]; and specialized methods including HUNTRESS [@HUNTRESS], SCIPhI [@SCIPhI], and Scelestial [@Scelestial].

However, the software landscape remains highly fragmented. Existing methods typically function as standalone binaries or scripts with inconsistent input/output formats, rendering comparative analysis difficult. Researchers attempting to utilize these tools face a laborious process of installation, data conversion, and script wrapping. Furthermore, integrating these tools into modern Python-based environments (e.g., alongside SCANPY [@SCANPY] or Biopython [@Biopython]) often requires custom development. `scphylo-tools` addresses these challenges by providing a single, cohesive Python API that wraps independently developed standalone applications, including Trisicell [@Trisicell-Boost], SCITE, PhISCS, and HUNTRESS.

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

`scphylo-tools` follows a modular architecture inspired by SCANPY. Modules mirror the analysis workflow: `io` (genotype and tree I/O, VCF export), `pp` (filtering of binary and read-count matrices, node collapsing, subsampling), `tl` (solvers, tree-comparison scores, consensus trees, partition functions, ancestral-state reconstruction), `pl` (tree rendering and genotype heatmaps), `ul` (utilities), and `datasets`. A `scphylo` console script exposes the same operations to shell users and workflow managers.

A single data model spans the package: genotype matrices are pandas `DataFrames` encoding presence, absence, and missing entries; read counts are AnnData objects; trees are NetworkX `DiGraph` objects. Over twenty solvers share the call shape `scp.tl.<method>(df_input, ...)` and return a conflict-free matrix, so switching methods is a one-line change; serialization, execution, parsing, and logging are internal. The ten tree-comparison metrics take the same two arguments, reducing a benchmark to a short loop.

To remove the main installation barrier, performance-critical C++ codebases (SCITE, ScisTree, MLTD, and the partition-function sampler) are vendored, compiled as Cython extensions, and shipped as prebuilt wheels. Non-redistributable tools stay optional, located by one resolution layer (explicit path, `SCPHYLO_TOOLS_DIR`, `scphylo.settings.tools_dir`, `PATH`) raising typed exceptions that name the missing artifact and how to supply it. Heavy optional dependencies (Gurobi, rpy2, graph-tool, Graphviz) load lazily.

Adding a method requires one module under `tl/solver` and no changes elsewhere, as with the existing wrappers. For validation, `datasets` retrieves published SCS datasets by identifier and the simulation engine generates ground-truth genotypes under configurable false-positive, false-negative, missing-value, and doublet noise.

# Research Impact Statement

`scphylo-tools` has enabled published research in three ways.

First, it provided the interface through which Trisicell and ScisTree were applied to build tumor phylogenies. These reconstructions underpinned the analysis of clonal sublines that uncovered heterogeneity-driven immunotherapy resistance in melanoma [@Gruen_2023], the comparison of genomic and epigenomic evolutionary trajectories in long-read sequencing of single cell-derived subclones [@Liu_2025], the stochastic modeling that quantified the non-genomic contribution to subclone evolution [@Hirsch_2025], and the recovery of metastatic migration histories from single-cell methylation sequencing [@Sgootr].

Second, it supported the identification of expressed mutation profiles in single cells, which revealed subclonal expansion patterns and the therapeutic impact of intratumor heterogeneity [@Trisicell-Boost].

Third, its uniform solver interface, simulation engine, and tree-comparison metrics were used to benchmark competing inference methods on common inputs during the evaluation of HUNTRESS [@HUNTRESS].

The package therefore serves as both the analysis layer for biological studies and the comparison harness for methods development. It is distributed on PyPI with source on GitHub and documentation on Read the Docs.

# AI Usage Disclosure

No generative AI tools were used in the development of the `scphylo-tools` software, the writing of this manuscript, or the preparation of supporting documentation.

# Acknowledgements

We acknowledge the contributions of the open-source community and the authors of the underlying algorithms whose work is wrapped within this toolkit. This work was supported by the Intramural Research Program of the National Institutes of Health, National Cancer Institute.

# References
