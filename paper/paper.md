---
title: "scphylo-tools: A comprehensive Python toolkit for single-cell tumor phylogenetic analysis"
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
date: 24 December 2025
bibliography: paper.bib
---

# Summary

`scphylo-tools` is a comprehensive Python library designed to unify single-cell tumor phylogeny inference methods. It addresses the lack of standardization by providing a cohesive interface for data processing, tree reconstruction, and visualization. By streamlining these tasks, `scphylo` empowers researchers with limited programming expertise to easily install and utilize complex inference methods. It accelerates algorithm development, ensures reproducible benchmarking, and democratizes access to advanced computational oncology tools.

# Introduction

Cancer is a dynamic evolutionary process driven by the acquisition of somatic mutations and the competitive selection of clonal populations. The rapid development of Single-Cell Sequencing (SCS) technologies has enabled the profiling of genomic alterations at cellular resolution, offering unprecedented insights into intratumor heterogeneity [@Navin_2011; @Wang_2014]. Elucidating the phylogenetic relationships between tumor cells is essential for understanding metastasis [@Leung_2017; @Roper_2020], drug resistance [@Kim_2014; @Gopalan_2021; @Gruen_2023], and clonal dynamics [@Laks_2019; @Liu_2025; @Hirsch_2025].

However, reconstructing the evolutionary history (often called phylogeny) of a tumor from SCS data presents unique computational challenges. These datasets are often plagued by high noise levels, including Allele Drop-Out (ADO), false positives, missing data, and doublet artifacts [@ReviewBinary; @RashidiMehrabadi2022]. Consequently, a diverse array of computational tools has been developed to address these challenges.

# Statement of need

Although the bioinformatics community has produced a plethora of algorithms to infer tumor trees, the software landscape remains highly fragmented. Existing methods typically function as standalone binaries or scripts with inconsistent input/output formats, rendering comparative analysis difficult. These include:

- **Stochastic and Probabilistic methods:** Tools like SCITE [@SCITE], infSCITE [@infSCITE], OncoNEM [@OncoNEM], SiFit [@SiFit], SiCloneFit [@SiCloneFit], and B-SCITE [@B-SCITE].
- **Combinatorial and Constraint-based methods:** Tools such as PhISCS [@PhISCS], PhISCS-BnB [@PhISCS-BnB], SPhyR [@SPhyR], ScisTree [@ScisTree], Trisicell [@Trisicell-Boost], Sgootr[@Sgootr], gpps [@gpps], and SASC [@SASC].
- **Specialized inference approaches:** Including HUNTRESS [@HUNTRESS], SCIPhI [@SCIPhI], TRaIT [@TRaIT], Phyolin [@Phyolin], scVILP [@scVILP], Scelestial [@Scelestial], and generative models like GRMT [@GRMT].

Consequently, researchers attempting to utilize these tools face a laborious process of installation, data conversion, and script wrapping. Furthermore, integrating these tools into modern Python-based environments (e.g., alongside SCANPY [@SCANPY] or Biopython [@Biopython]) often requires custom development. `scphylo` addresses these challenges by wrapping diverse state-of-the-art algorithms into a single, cohesive Python API.

# Features and Functionality

The `scphylo` package is modular, catering to different stages of the phylogenetic analysis workflow:

## 1. Data Management (`scphylo.io`)

`scphylo` streamlines data management through its `io` module, providing a unified interface for reading and writing genotype and read-count matrices. This ensures seamless integration between upstream variant callers and downstream phylogenetic analysis, directly addressing the format inconsistencies common in the field.

## 2. Data Preprocessing (`scphylo.pp`)

To mitigate the noise and sparsity characteristic of SCS data, the `pp` module provides essential preprocessing utilities:

- **Filtering:** Functions to remove artifacts and retain informative features based on mutation prevalence.
- **Matrix Optimization:** The `bifiltering` algorithm identifies the maximally informed submatrix, minimizing the impact of missing data on tree reconstruction.
- **Data Aggregation:** The `consensus` utility enables cell merging to generate high-confidence genotypes or to cluster similar cells, significantly improving the signal-to-noise ratio.

## 3. Phylogenetic Inference (`scphylo.tl`)

The core of `scphylo` lies in its unified solver interface. It wraps a wide array of the inference algorithms mentioned above (e.g., SCITE, PhISCS, HUNTRESS, ScisTree), allowing users to execute complex inference tasks using standard Python syntax. This abstraction layer handles the formatting of input files, execution of the underlying binary/script, and parsing of the resulting trees back into a NetworkX-compatible format. This is critical for users wishing to compare different inference methods on the same data.

## 4. Tree Visualization (`scphylo.pl`)

Interpreting phylogenetic trees requires high-quality visualization. `scphylo` provides plotting functions that leverage Matplotlib and Graphviz to render annotated trees and sorted genotype heatmaps. These visualizations help researchers identify clonal subpopulations, visualize the acquisition of mutations, and characterize the evolutionary trajectories of distinct cell lineages.

## 5. Interactive Datasets (`scphylo.datasets`)

To foster reproducibility and ease of benchmarking, `scphylo` includes a dedicated module providing direct access to a curated collection of published high-impact SCS datasets. Users can load genotype matrices from landmark cancer evolution studies with a single command. The available datasets cover a wide range of cancer types and sequencing protocols, including:

- **Leukemia:** Acute Lymphoblastic Leukemia (ALL) [@Gawad_2014], Acute Myeloid Leukemia (AML) [@Morita_2020], and JAK2-negative Myeloproliferative Neoplasms [@Hou_2012].
- **Solid Tumors:** Kidney tumor (single-cell exome) [@Xu_2012], Bladder cancer [@Li_2012], Breast cancer [@Wang_2014], and Melanoma [@Wolf_2019].
- **Complex Evolution:** Metastatic Colorectal Cancer [@Leung_2017], Non-hereditary Colorectal Cancer [@Wu_2016], High-grade Serous Ovarian Cancer [@McPherson_2016; @McPherson_2019], and Oligodendroglioma [@Tirosh_2016].
- **Protocol-Specific Data:** Datasets generated via SNES [@Leung_2015] and approaches for identifying leukemic stem cells [@Velten_2021].

## 6. Tree Evaluation and Simulation (`scphylo.ul`)

To facilitate rigorous benchmarking, `scphylo` implements a comprehensive suite of metrics to quantitatively assess the accuracy of inferred trees relative to a ground truth. Supported metrics include:

- **Lineage Accuracies:** Assessments of how closely predicted lineage relationships match the ground truth, including Ancestor-Descendant, Different-Lineage, Robinson-Foulds, and Genotype-Similarity.
- **Tree Distance Measures:** Implementation of specialized distances such as the Multi-Labeled Tree Dissimilarity (MLTD) [@MLTD], CASet and DISC distances [@CASet_DISC], MP3 similarity score [@MP3], and Bourque distances [@Bourque].
- **Subclonal Reliability:** The Trisicell-PartF score [@Trisicell-PartF], which utilizes a partition function algorithm to evaluate the reliability of inferred subclonal structures.

To complement these benchmarking tools, `scphylo` features a flexible simulation engine for generating synthetic tumor phylogenies with known ground truths. This is critical for validating new algorithms and profiling performance across different noise regimes. Users can configure:

- **Topology:** Parameters for defining tree structure, including tree size, branching factors, and mutation attachment strategies.
- **Noise:** Injection of realistic technical errors with customizable rates for false positives, false negatives (ADO), and missing data.
- **Artifacts:** Simulation of doublet cells to mimic common confounders found in SCS data.

# Acknowledgements

We acknowledge the contributions of the open-source community and the authors of the underlying algorithms mentioned throughout this paper, whose work is wrapped or utilized within this toolkit.

# References
