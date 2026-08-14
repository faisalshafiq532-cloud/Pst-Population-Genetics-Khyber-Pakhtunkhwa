# ============================================================
# Isolation-by-Distance Analysis Using a Mantel Test
#
# Study: Puccinia striiformis f. sp. tritici populations from
#        Khyber Pakhtunkhwa, Pakistan
#
# Inputs:
#   SSR_PK-data.xlsx
#   Locations_Sampling.xlsx
#
# Analyses:
#   1. Bruvo genetic distance
#   2. Geographic distance based on the Haversine method
#   3. Mantel test for isolation by distance
#
# Output:
#   IBD_Mantel/
#     IBD_pairwise_distances.csv
#     IBD_Mantel_results.txt
#     IBD_Mantel_plot.tiff
#
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(readxl)
library(adegenet)
library(poppr)
library(geosphere)
library(vegan)
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
# This conversion follows the genotype coding used in the
# original analysis workflow.
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
# 4. Create genind object
# ------------------------------------------------------------

gen <- df2genind(
  ssr_genotypes,
  sep = "/",
  ploidy = 2,
  ind.names = paste0(
    "PK_",
    seq_len(nrow(ssr_data))
  ),
  pop = ssr_data$Regions,
  type = "codom"
)


# ------------------------------------------------------------
# 5. Define SSR repeat lengths
#
# Values correspond to the 17 SSR loci used in the analysis.
# ------------------------------------------------------------

repeat_lengths <- c(
  nrjn12 = 3,
  nrjn8  = 3,
  nrjn13 = 2,
  nrjn3  = 2,
  nrjn11 = 2,
  nrjo27 = 2,
  nrjn6  = 3,
  nrjo21 = 3,
  nrjn10 = 3,
  nrjo18 = 3,
  nwu6   = 2,
  nrjo20 = 3,
  nrjn2  = 2,
  nrjn4  = 2,
  nrjn9  = 2,
  nrjn5  = 2,
  nwu12  = 3
)

stopifnot(
  identical(
    names(repeat_lengths),
    locNames(gen)
  )
)


# ============================================================
# PART I: GENETIC DISTANCE
# ============================================================


# ------------------------------------------------------------
# 6. Calculate Bruvo genetic distance
# ------------------------------------------------------------

bruvo_dist <- poppr::bruvo.dist(
  gen,
  replen = repeat_lengths,
  add = TRUE
)

stopifnot(
  length(labels(bruvo_dist)) == nInd(gen)
)


# ============================================================
# PART II: GEOGRAPHIC DISTANCE
# ============================================================


# ------------------------------------------------------------
# 7. Import sampling coordinates
# ------------------------------------------------------------

coords <- read_excel(
  "Locations_Sampling.xlsx"
)


# ------------------------------------------------------------
# 8. Check coordinate dataset size
#
# The original analysis used coordinate records in the same
# isolate order as the genind object.
# ------------------------------------------------------------

if (nrow(coords) != nInd(gen)) {
  
  stop(
    "The number of coordinate records does not match ",
    "the number of isolates."
  )
}


# ------------------------------------------------------------
# 9. Assign isolate identifiers
# ------------------------------------------------------------

coords$Isolate <- indNames(gen)


# ------------------------------------------------------------
# 10. Extract geographic coordinates
# ------------------------------------------------------------

geo_coords <- coords[
  ,
  c(
    "Longitude (E)",
    "Latitude (N)"
  )
]

geo_coords <- as.matrix(
  geo_coords
)


# ------------------------------------------------------------
# 11. Calculate pairwise geographic distances
#
# Haversine distances are calculated in metres and converted
# to kilometres.
# ------------------------------------------------------------

geo_dist_matrix <- geosphere::distm(
  geo_coords,
  fun = geosphere::distHaversine
)

geo_dist_matrix <- geo_dist_matrix / 1000

rownames(geo_dist_matrix) <- coords$Isolate
colnames(geo_dist_matrix) <- coords$Isolate

geo_dist <- as.dist(
  geo_dist_matrix
)


# ------------------------------------------------------------
# 12. Verify distance matrix compatibility
# ------------------------------------------------------------

if (
  !identical(
    labels(bruvo_dist),
    labels(geo_dist)
  )
) {
  
  stop(
    "Genetic and geographic distance matrices do not have ",
    "matching isolate labels."
  )
}


# ============================================================
# PART III: MANTEL TEST
# ============================================================


# ------------------------------------------------------------
# 13. Test isolation by distance
#
# Pearson Mantel correlation with 9,999 permutations.
# ------------------------------------------------------------

set.seed(20260809)

mantel_ibd <- vegan::mantel(
  bruvo_dist,
  geo_dist,
  method = "pearson",
  permutations = 9999
)


# ============================================================
# PART IV: PREPARE DATA FOR VISUALIZATION
# ============================================================


# ------------------------------------------------------------
# 14. Extract pairwise genetic and geographic distances
# ------------------------------------------------------------

ibd_df <- data.frame(
  Genetic_distance = as.vector(
    bruvo_dist
  ),
  Geographic_distance = as.vector(
    geo_dist
  )
)


# ============================================================
# PART V: CREATE MANUSCRIPT FIGURE
# ============================================================


# ------------------------------------------------------------
# 15. Define Mantel test statistics for annotation
# ------------------------------------------------------------

mantel_r <- round(
  unname(mantel_ibd$statistic),
  3
)

mantel_p <- mantel_ibd$signif

n_pairs <- nrow(
  ibd_df
)


# ------------------------------------------------------------
# 16. Create IBD scatter plot
# ------------------------------------------------------------

ibd_plot <- ggplot(
  ibd_df,
  aes(
    x = Geographic_distance,
    y = Genetic_distance
  )
) +
  
  geom_point(
    color = "#2C7FB8",
    alpha = 0.20,
    size = 1.2
  ) +
  
  geom_smooth(
    method = "lm",
    color = "#D95F02",
    fill = "#D95F02",
    alpha = 0.20,
    linewidth = 1
  ) +
  
  annotate(
    "text",
    x = max(
      ibd_df$Geographic_distance,
      na.rm = TRUE
    ) * 0.65,
    y = max(
      ibd_df$Genetic_distance,
      na.rm = TRUE
    ) * 0.15,
    label = paste0(
      "Mantel r = ",
      sprintf("%.3f", mantel_r),
      "\nP ",
      ifelse(
        mantel_p < 0.001,
        "< 0.001",
        paste0(
          "= ",
          sprintf("%.3f", mantel_p)
        )
      ),
      "\nn = ",
      format(
        n_pairs,
        big.mark = ","
      )
    ),
    size = 5,
    hjust = 0
  ) +
  
  labs(
    x = "Geographic distance (km)",
    y = "Genetic distance (Bruvo distance)"
  ) +
  
  theme_classic(
    base_size = 16
  ) +
  
  theme(
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_text(
      color = "black"
    ),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    axis.line = element_blank()
  )


# ============================================================
# PART VI: EXPORT RESULTS
# ============================================================


# ------------------------------------------------------------
# 17. Create output directory
# ------------------------------------------------------------

output_dir <- "IBD_Mantel"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# ------------------------------------------------------------
# 18. Export pairwise distance data
# ------------------------------------------------------------

write.csv(
  ibd_df,
  file = file.path(
    output_dir,
    "IBD_pairwise_distances.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 19. Export Mantel test results
# ------------------------------------------------------------

capture.output(
  mantel_ibd,
  file = file.path(
    output_dir,
    "IBD_Mantel_results.txt"
  )
)


# ------------------------------------------------------------
# 20. Export publication-quality TIFF figure
# ------------------------------------------------------------

tiff(
  filename = file.path(
    output_dir,
    "IBD_Mantel_plot.tiff"
  ),
  width = 7.5,
  height = 5,
  units = "in",
  res = 600,
  compression = "lzw"
)

print(
  ibd_plot
)

dev.off()


# ------------------------------------------------------------
# 21. Report results
# ------------------------------------------------------------

cat(
  "\nIsolation-by-distance analysis complete.\n\n"
)

cat(
  "Number of isolates:",
  nInd(gen),
  "\n"
)

cat(
  "Number of pairwise comparisons:",
  n_pairs,
  "\n"
)

cat(
  "Mantel r:",
  sprintf("%.3f", mantel_r),
  "\n"
)

cat(
  "Mantel P-value:",
  mantel_p,
  "\n"
)

print(
  mantel_ibd
)