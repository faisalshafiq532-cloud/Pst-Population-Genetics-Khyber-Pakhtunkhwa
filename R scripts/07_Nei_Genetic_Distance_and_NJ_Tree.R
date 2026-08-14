# ============================================================
# Population-Level Nei Genetic Distance and Neighbor-Joining Tree
#
# Study: Puccinia striiformis f. sp. tritici populations from
#        Khyber Pakhtunkhwa, Pakistan
#
# Input:
#   SSR_PK-data.xlsx
#   Sheet: PK
#
# Analyses:
#   1. Population-level Nei genetic distance
#   2. Neighbor-Joining tree based on Nei genetic distance
#
# Outputs:
#   Nei_NJ_Tree/
#     Nei_Genetic_Distance_Matrix.csv
#     Nei_NJ_Tree_Population.tiff
#
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(readxl)
library(adegenet)
library(poppr)
library(ape)
library(ggtree)
library(ggplot2)


# ------------------------------------------------------------
# 2. Import SSR genotype data
# ------------------------------------------------------------

input_file <- "SSR_PK-data.xlsx"

ssr_data <- read_excel(
  input_file,
  sheet = "PK"
)

# The 17 SSR loci are located in columns 4-20.
locus_cols <- names(ssr_data)[4:20]

ssr_genotypes <- ssr_data[locus_cols]


# ------------------------------------------------------------
# 3. Convert SSR genotypes to allele-pair format
#
# This follows the genotype coding used in the original analysis.
# ------------------------------------------------------------

for (locus in locus_cols) {
  
  genotype <- as.character(
    ssr_genotypes[[locus]]
  )
  
  allele1 <- substr(
    genotype,
    1,
    3
  )
  
  allele2 <- substr(
    genotype,
    4,
    6
  )
  
  ssr_genotypes[[locus]] <- paste(
    allele1,
    allele2,
    sep = "/"
  )
}


# ------------------------------------------------------------
# 4. Define manuscript population order
# ------------------------------------------------------------

district_order <- c(
  "Peshawar",
  "Kohat",
  "Charsadda",
  "Mardan",
  "Swabi",
  "Manshera",
  "Buner"
)


# ------------------------------------------------------------
# 5. Create genind object
# ------------------------------------------------------------

gen <- df2genind(
  ssr_genotypes,
  sep = "/",
  ploidy = 2,
  ind.names = paste0(
    "PK_",
    seq_len(nrow(ssr_data))
  ),
  pop = factor(
    ssr_data$Regions,
    levels = district_order
  ),
  type = "codom"
)


# ============================================================
# PART I: POPULATION-LEVEL NEI GENETIC DISTANCE
# ============================================================


# ------------------------------------------------------------
# 6. Aggregate individual genotypes into populations
# ------------------------------------------------------------

genpop_obj <- genind2genpop(
  gen
)


# ------------------------------------------------------------
# 7. Calculate Nei genetic distance among populations
# ------------------------------------------------------------

nei_pop_dist <- poppr::nei.dist(
  genpop_obj
)

# Convert distance object to matrix.
nei_pop_matrix <- as.matrix(
  nei_pop_dist
)

# Arrange populations in manuscript order.
population_order <- district_order[
  district_order %in% rownames(nei_pop_matrix)
]

nei_pop_matrix <- nei_pop_matrix[
  population_order,
  population_order,
  drop = FALSE
]

# Recreate distance object after ordering.
nei_pop_dist <- as.dist(
  nei_pop_matrix
)


# ============================================================
# PART II: NEIGHBOR-JOINING TREE
# ============================================================


# ------------------------------------------------------------
# 8. Construct Neighbor-Joining tree
# ------------------------------------------------------------

nei_nj_tree <- ape::nj(
  nei_pop_dist
)


# ------------------------------------------------------------
# 9. Create visualization copy of the NJ tree
#
# Negative branch lengths, if present, are set to zero only for
# visualization. The original tree remains unchanged.
# ------------------------------------------------------------

nei_nj_tree_plot <- nei_nj_tree

if (any(nei_nj_tree_plot$edge.length < 0)) {
  
  nei_nj_tree_plot$edge.length[
    nei_nj_tree_plot$edge.length < 0
  ] <- 0
}


# ============================================================
# PART III: TREE VISUALIZATION
# ============================================================


# ------------------------------------------------------------
# 10. Define population colours
# ------------------------------------------------------------

district_colors <- c(
  "Peshawar"  = "#1B9E77",
  "Kohat"     = "#D95F02",
  "Charsadda" = "#7570B3",
  "Mardan"    = "#E7298A",
  "Swabi"     = "#66A61E",
  "Manshera"  = "#E6AB02",
  "Buner"     = "#A6761D"
)


# ------------------------------------------------------------
# 11. Create circular Neighbor-Joining tree
# ------------------------------------------------------------

nei_tree_plot <- ggtree(
  nei_nj_tree_plot,
  layout = "circular",
  linewidth = 1.2
) +
  
  geom_tippoint(
    aes(color = label),
    size = 4
  ) +
  
  geom_tiplab(
    aes(color = label),
    size = 5,
    fontface = "bold",
    offset = 0.045,
    align = FALSE
  ) +
  
  scale_color_manual(
    values = district_colors
  ) +
  
  theme(
    legend.position = "none"
  )


# ============================================================
# PART IV: EXPORT RESULTS
# ============================================================


# ------------------------------------------------------------
# 12. Create output directory
# ------------------------------------------------------------

output_dir <- "Nei_NJ_Tree"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# ------------------------------------------------------------
# 13. Export Nei genetic distance matrix
# ------------------------------------------------------------

write.csv(
  nei_pop_matrix,
  file = file.path(
    output_dir,
    "Nei_Genetic_Distance_Matrix.csv"
  ),
  row.names = TRUE
)


# ------------------------------------------------------------
# 14. Export Neighbor-Joining tree in Newick format
# ------------------------------------------------------------

write.tree(
  nei_nj_tree,
  file = file.path(
    output_dir,
    "Nei_NJ_Tree.newick"
  )
)


# ------------------------------------------------------------
# 15. Export publication-quality TIFF figure
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "Nei_NJ_Tree_Population.tiff"
  ),
  plot = nei_tree_plot,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


# ------------------------------------------------------------
# 16. Report key results
# ------------------------------------------------------------

cat(
  "\nNei genetic distance and NJ tree analysis complete.\n\n"
)

cat(
  "Number of populations:",
  nPop(genpop_obj),
  "\n"
)

cat(
  "Population names:",
  paste(
    population_order,
    collapse = ", "
  ),
  "\n"
)

cat(
  "\nMinimum Nei genetic distance:",
  round(
    min(
      nei_pop_matrix[
        upper.tri(nei_pop_matrix)
      ]
    ),
    4
  ),
  "\n"
)

cat(
  "Maximum Nei genetic distance:",
  round(
    max(
      nei_pop_matrix[
        upper.tri(nei_pop_matrix)
      ]
    ),
    4
  ),
  "\n"
)

cat(
  "\nNJ tree exported to:",
  output_dir,
  "\n"
)

print(
  round(
    nei_pop_matrix,
    4
  )
)

print(
  nei_tree_plot
)