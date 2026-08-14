# Analysis Scripts

This directory contains the cleaned R scripts used for the population genetic, spatial, and connectivity analyses described in the associated manuscript.

The scripts were organized from the original analysis workflow for public release.

## Scripts

### 01. SSR marker quality and diversity

`01_SSR_Marker_Quality_and_Table1.R`

Assessment of SSR marker quality and diversity and generation of the corresponding manuscript table.

### 02. Population genetic and genotypic diversity

`02_Population_Genetic_and_Genotypic_Diversity.R`

Analysis of population genetic diversity and genotypic diversity across the sampled populations.

### 03. Multilocus genotype distribution

`03_MLG_Distribution_Table2.R`

Identification and distribution of multilocus genotypes (MLGs) and generation of the corresponding manuscript table.

### 04. MLG distribution visualization

`04_MLG_Distribution_Visualization.R`

Visualization of the distribution and frequency of multilocus genotypes across populations.

### 05. AMOVA and pairwise FST

`05_AMOVA_and_Pairwise_FST.R`

Analysis of molecular variance (AMOVA) and pairwise Weir-Cockerham FST among populations.

### 06. Isolation by distance

`06_Isolation_by_Distance_Mantel.R`

Isolation-by-distance analysis using Bruvo genetic distance, geographic distance, and Mantel tests.

### 07. Nei's genetic distance and neighbor-joining analysis

`07_Nei_Genetic_Distance_NJ_Tree.R`

Calculation of Nei's genetic distance among populations and construction of a neighbor-joining tree.

### 08. Spatial MLG distribution

`08_Spatial_MLG_Distribution.R`

Spatial analysis and visualization of multilocus genotype distribution.

### 09. Spatial visualization of estimated historical gene flow

`09_Spatial_Visualization_of_Estimated_Gene_Flow.R`

Spatial visualization of estimated historical gene flow among populations.

### 10. Principal coordinates analysis

`10_PCoA_Population_Structure.R`

Principal coordinates analysis (PCoA) based on genetic distance.

### 11. Discriminant analysis of principal components

`11_DAPC_Population_Structure.R`

Discriminant analysis of principal components (DAPC) to examine population structure.

### 12. STRUCTURE visualization

`12_STRUCTURE_Plotting.R`

Visualization of population structure inferred using STRUCTURE.

STRUCTURE analyses and alignment of replicate runs were conducted using external software. The R script uses the corresponding prepared output files to generate the population structure visualization.

### 13. Atmospheric connectivity analysis

`13_Atmospheric_Connectivity_Analysis.R`

Analysis of atmospheric connectivity using external trajectory and atmospheric data sources.

The underlying external datasets and downloaded trajectory data are not redistributed in this repository. The script is provided to document the computational workflow used for the analysis.

## Input data

The primary SSR genotype dataset is located at:

`../data/SSR_PK-data.xlsx`

Files required for the STRUCTURE visualization workflow are located at:

`../data/STRUCTURE/`

Users should ensure that the required input files are available before running the corresponding scripts.

## External software

Selected analyses relied on external software in addition to R, including STRUCTURE and CLUMPP.

The atmospheric connectivity workflow also relies on external trajectory and atmospheric datasets that are not redistributed in this repository.

## Reproducibility

The scripts were cleaned from the original working analysis files for public release while retaining the computational procedures used for the reported analyses.

Some analyses depend on objects, metadata, or results generated during preceding analytical steps. Users should therefore review the comments and input requirements at the beginning of each script before execution.
