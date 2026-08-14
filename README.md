# Population Genetic Analysis of *Puccinia striiformis* f. sp. *tritici* from Khyber Pakhtunkhwa, Pakistan

This repository contains the R scripts and supporting data used for population genetic analyses of 144 single-urediniospore-derived isolates of *Puccinia striiformis* f. sp. *tritici* collected from seven districts of Khyber Pakhtunkhwa, Pakistan and genotyped using 17 SSR markers.

The repository accompanies the manuscript describing genetic diversity, multilocus genotype distribution, population differentiation, spatial structure, and connectivity of the *P. striiformis* f. sp. *tritici* population.

## Repository structure

```text
├── R scripts/
│   ├── 01_SSR_Marker_Quality_and_Table1.R
│   ├── 02_Population_Genetic_and_Genotypic_Diversity.R
│   ├── 03_MLG_Distribution_Table2.R
│   ├── 04_MLG_Distribution_Figure.R
│   ├── 05_AMOVA_and_Pairwise_FST.R
│   ├── 06_Isolation_by_Distance_Mantel.R
│   ├── 07_Nei_Genetic_Distance_NJ_Tree.R
│   ├── 08_Spatial_MLG_Distribution.R
│   ├── 09_Spatial_Visualization_of_Estimated_Gene_Flow.R
│   ├── 10_PCoA_Population_Structure.R
│   ├── 11_DAPC_Population_Structure.R
│   ├── 12_STRUCTURE_Plotting.R
│   └── 13_Atmospheric_Connectivity_Analysis.R
│
└── data/
    ├── SSR_PK-data.xlsx
    └── STRUCTURE/
```

## Data

The file `data/SSR_PK-data.xlsx` contains the SSR genotype dataset used for the population genetic analyses.

The `data/STRUCTURE/` directory contains files required for the visualization of population structure results.

Supplementary tables and figures associated with the manuscript are provided through the journal as supplementary material and are not duplicated in this repository.

## Analyses

The scripts include analyses of:

1. SSR marker quality and diversity
2. Population genetic and genotypic diversity
3. Multilocus genotype (MLG) distribution
4. MLG distribution visualization
5. AMOVA and pairwise Weir-Cockerham FST
6. Isolation by distance using Mantel tests
7. Nei's genetic distance and neighbor-joining analysis
8. Spatial MLG distribution
9. Spatial visualization of estimated historical gene flow
10. Principal coordinates analysis (PCoA)
11. Discriminant analysis of principal components (DAPC)
12. Visualization of population structure inferred using STRUCTURE
13. Atmospheric connectivity analysis

## Software and R packages

Analyses were conducted in R using packages including:

- adegenet
- ape
- dplyr
- geosphere
- ggplot2
- ggrepel
- ggtree
- hierfstat
- openxlsx
- pegas
- poppr
- readxl
- vegan
- writexl

Additional external software was used for selected analyses, including STRUCTURE and CLUMPP. The corresponding R script in this repository was used for visualization of population structure results.

## Atmospheric connectivity analysis

The atmospheric connectivity analysis was conducted using external trajectory and atmospheric data sources. Because the underlying external datasets and downloaded trajectory data are not redistributed in this repository, the repository provides the cleaned analysis script documenting the computational workflow.

## Reproducibility

The scripts were cleaned and organized from the original analysis workflow for public release. Users should review the required packages and input files before running individual scripts.

Some analyses depend on objects or results generated during preceding analytical steps. Scripts should therefore be run with attention to their stated input requirements and analysis sequence.

## Citation

If you use the scripts or data from this repository, please cite the associated manuscript:

[Full citation will be added after publication.]
