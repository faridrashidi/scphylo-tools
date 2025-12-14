---
title: 'scphylo: A comprehensive Python toolkit for single-cell tumor phylogenetic analysis'
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
  - name: Cancer Data Science Laboratory, Center for Cancer Research, National Cancer Institute, National Institutes of Health, Bethesda, MD, USA
    index: 1
  - name: Laboratory of Human Carcinogenesis, Center for Cancer Research, National Cancer Institute, National Institutes of Health, Bethesda, MD, USA
    index: 2
date: 13 December 2025
bibliography: paper.bib
---

# Summary

Cancer is an evolutionary process characterized by the accumulation of somatic mutations and the competitive expansion of clonal populations. The advent of Single-Cell Sequencing (SCS) has revolutionized our ability to study this process by providing high-resolution data on the genomic profiles of individual tumor cells. However, inferring the evolutionary history (phylogeny) of a tumor from SCS data is computationally challenging due to significant technical noise, including Allele Drop-Out (ADO), false positives, and doublets, as well as the complexity of missing data.

`scphylo` is a Python library designed to facilitate the analysis, simulation, and visualization of single-cell tumor phylogenies. It provides a unified interface for manipulating single-cell genotypes, simulating realistic tumor data with various noise profiles, and inferring phylogenetic trees using state-of-the-art algorithms. By bridging the gap between data generation, processing, and inference, `scphylo` aims to streamline computational oncology workflows and enhance the reproducibility of cancer evolution studies.

# Statement of need

Understanding the evolutionary history of tumors is critical for identifying drivers of metastasis and mechanisms of drug resistance. While numerous algorithms have been developed to infer tumor phylogenies from noisy SCS data—such as SCITE [@jahn2016accurate], PhISCS [@malikic2019phiscs], and HUNTRESS [@govek2020huntress]—the ecosystem for these tools remains fragmented.

Many existing methods are distributed as standalone binaries or scripts in varying languages (C++, Java, R) with inconsistent input/output formats. This fragmentation creates significant barriers for researchers who need to:
1.  **Benchmark algorithms:** Comparing the performance of different inference methods on a standardized dataset is currently a laborious manual process.
2.  **Pipeline integration:** Integrating phylogeny inference into larger Python-based bioinformatics pipelines (e.g., alongside `Scanpy` or `Biopython`) requires writing custom wrappers for each tool.
3.  **Visualization:** High-quality, programmatic visualization of tumor trees with mutation annotations is often lacking in solver-specific packages.

`scphylo` addresses these challenges by providing a comprehensive "gym" for single-cell phylogeny. It offers a standardized Python API that wraps multiple inference algorithms, allowing users to run and compare methods like SCITE, PhISCS, and HUNTRESS with a few lines of code. Furthermore, it includes robust utilities for generating synthetic datasets with ground-truth trees, which are essential for validating new methods. By centralizing these capabilities, `scphylo` reduces the friction in computational oncology research and democratizes access to advanced phylogenetic analysis tools.

# Features and Functionality

The `scphylo` package is structured into several modular components, ensuring flexibility and ease of use:

## 1. Data Input/Output and Preprocessing (`scphylo.io`)
Real-world SCS data comes in various formats. `scphylo` supports reading and writing genotypes in standard formats (e.g., FASTA, CSV, GML) and converts them into an internal representation optimized for analysis. It handles the parsing of mutation matrices where rows represent cells and columns represent loci, automatically managing missing data and ternary states (0: homozygous reference, 1: heterozygous, 2: homozygous variant).

## 2. Phylogenetic Inference (`scphylo.tl`)
A core feature of the library is its solver interface. `scphylo` provides wrappers for several popular phylogenetic inference tools. Users can infer trees using:
* **PhISCS:** A combinatorial approach for perfect phylogeny using constraint satisfaction.
* **HUNTRESS:** A fast reconstruction algorithm designed to handle high ADO rates.
* **SCITE:** A stochastic search algorithm (MCMC) that handles noisy data.

This unified interface allows researchers to run ensemble methods or compare the topology of trees inferred by different algorithms on the same dataset.

## 3. Data Simulation (`scphylo.ul`)
To rigorously test inference algorithms, researchers need ground-truth data. `scphylo` includes a sophisticated simulation engine that generates synthetic tumor phylogenies and genotype matrices. Users can control:
* **Tree Topology:** Number of cells, number of mutations, and branching factors.
* **Noise Models:** Customizable rates for Allele Drop-Out (false negatives), false positives, and missing data entries.
* **Doublets:** Simulation of doublet artifacts common in single-cell protocols.

## 4. Visualization (`scphylo.pl`)
Interpreting the output of phylogenetic inference requires clear visualization. `scphylo` leverages `networkx` and `matplotlib` to produce publication-quality plots of phylogenetic trees. Users can visualize the clonal evolution, annotate nodes with specific mutations, and generate heatmaps of the genotype matrices sorted by phylogenetic placement, making it easier to identify exclusive sub-clones.

## 5. Metrics and Benchmarking
The library implements standard metrics for comparing phylogenetic trees, such as the Ancestor-Descendant accuracy and multi-state pairwise distances. These metrics enable quantitative benchmarking of reconstructed trees against the ground truth.

# Availability

`scphylo` is open-source software released under the BSD-3-Clause license. It is available on the Python Package Index (PyPI) and can be installed via `pip`. The source code, documentation, and issue tracker are hosted on GitHub.

# Acknowledgements

We acknowledge the contributions of the open-source community and the authors of the underlying algorithms (SCITE, PhISCS, HUNTRESS) whose work is wrapped or utilized within this toolkit.

# References
