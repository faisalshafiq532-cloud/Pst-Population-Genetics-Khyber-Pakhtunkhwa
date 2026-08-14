# ============================================================
# Principal Coordinates Analysis (PCoA) Based on Bruvo Distance
#
# Study: Puccinia striiformis f. sp. tritici populations from
#        Khyber Pakhtunkhwa, Pakistan
#
# Required input objects:
#   bruvo_dist       : Bruvo genetic distance matrix
#   mlg_geo_correct  : isolate metadata containing Isolate,
#                      District, MLG, and MLG_group
#
# Analyses:
#   1. Perform PCoA on Bruvo genetic distances
#   2. Extract coordinates for the first two PCoA axes
#   3. Calculate the variance explained by each axis
#   4. Merge PCoA coordinates with population metadata
#   5. Calculate population centroids
#   6. Visualize isolates and population centroids
#
# Outputs:
#   PCoA_coordinates_population.csv
#   PCoA_population_structure.tiff
#
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(ape)
library(dplyr)
library(ggplot2)
library(ggrepel)


# ============================================================
# PART I: VALIDATE INPUT OBJECTS
# ============================================================


# ------------------------------------------------------------
# 2. Check required objects
# ------------------------------------------------------------

required_objects <- c(
  "bruvo_dist",
  "mlg_geo_correct"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1)
  )
]

if (length(missing_objects) > 0) {
  
  stop(
    "The following required objects are missing: ",
    paste(missing_objects, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 3. Check required metadata columns
# ------------------------------------------------------------

required_columns <- c(
  "Isolate",
  "District",
  "MLG",
  "MLG_group"
)

missing_columns <- required_columns[
  !required_columns %in% names(mlg_geo_correct)
]

if (length(missing_columns) > 0) {
  
  stop(
    "The following columns are missing from mlg_geo_correct: ",
    paste(missing_columns, collapse = ", ")
  )
}


# ============================================================
# PART II: PRINCIPAL COORDINATES ANALYSIS
# ============================================================


# ------------------------------------------------------------
# 4. Perform PCoA using Bruvo genetic distance
# ------------------------------------------------------------

pcoa_result <- ape::pcoa(
  bruvo_dist
)


# ------------------------------------------------------------
# 5. Extract coordinates for the first two PCoA axes
# ------------------------------------------------------------

pcoa_coords <- as.data.frame(
  pcoa_result$vectors[, 1:2]
)

colnames(pcoa_coords) <- c(
  "PCoA1",
  "PCoA2"
)

pcoa_coords$Isolate <- rownames(
  pcoa_coords
)


# ------------------------------------------------------------
# 6. Calculate variance explained by the first two axes
# ------------------------------------------------------------

pcoa_variance <- pcoa_result$values$Relative_eig[
  1:2
] * 100


# ============================================================
# PART III: MERGE PCOA COORDINATES WITH METADATA
# ============================================================


# ------------------------------------------------------------
# 7. Merge PCoA coordinates with isolate metadata
# ------------------------------------------------------------

pcoa_df <- pcoa_coords %>%
  inner_join(
    mlg_geo_correct %>%
      select(
        Isolate,
        District,
        MLG,
        MLG_group
      ),
    by = "Isolate"
  )


# ------------------------------------------------------------
# 8. Verify that all isolates were retained
# ------------------------------------------------------------

if (nrow(pcoa_df) != nrow(pcoa_coords)) {
  
  warning(
    "Not all PCoA isolates were matched to metadata."
  )
}


# ============================================================
# PART IV: DEFINE POPULATION COLOURS
# ============================================================


district_colors <- c(
  "Peshawar"  = "#1B9E77",
  "Kohat"     = "#D95F02",
  "Charsadda" = "#7570B3",
  "Mardan"    = "#E7298A",
  "Swabi"     = "#66A61E",
  "Manshera"  = "#E6AB02",
  "Buner"     = "#A6761D"
)


# ============================================================
# PART V: CALCULATE POPULATION CENTROIDS
# ============================================================


pcoa_centroids <- pcoa_df %>%
  group_by(
    District
  ) %>%
  summarise(
    PCoA1 = mean(PCoA1),
    PCoA2 = mean(PCoA2),
    n = n(),
    .groups = "drop"
  )


# ============================================================
# PART VI: CREATE PCOA PLOT
# ============================================================


pcoa_plot <- ggplot(
  pcoa_df,
  aes(
    x = PCoA1,
    y = PCoA2,
    color = District
  )
) +
  
  # Reference lines
  geom_hline(
    yintercept = 0,
    color = "grey85",
    linewidth = 0.3
  ) +
  
  geom_vline(
    xintercept = 0,
    color = "grey85",
    linewidth = 0.3
  ) +
  
  # Individual isolates
  geom_point(
    size = 3,
    alpha = 0.75
  ) +
  
  # Population centroids
  geom_point(
    data = pcoa_centroids,
    aes(
      x = PCoA1,
      y = PCoA2,
      color = District
    ),
    shape = 4,
    size = 6,
    stroke = 1.5
  ) +
  
  # Population labels
  geom_text_repel(
    data = pcoa_centroids,
    aes(
      x = PCoA1,
      y = PCoA2,
      label = District,
      color = District
    ),
    size = 4,
    fontface = "bold",
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  
  # Population colours
  scale_color_manual(
    values = district_colors
  ) +
  
  # Axis labels
  labs(
    x = paste0(
      "PCoA1 (",
      round(pcoa_variance[1], 2),
      "%)"
    ),
    y = paste0(
      "PCoA2 (",
      round(pcoa_variance[2], 2),
      "%)"
    ),
    color = "Population"
  ) +
  
  theme_classic(
    base_size = 15
  ) +
  
  theme(
    legend.position = "right",
    
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    
    legend.text = element_text(
      size = 11
    ),
    
    axis.title = element_text(
      size = 15,
      face = "bold"
    ),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    axis.ticks.length = grid::unit(
      0.15,
      "cm"
    )
  )


# ============================================================
# PART VII: EXPORT RESULTS
# ============================================================


# ------------------------------------------------------------
# 9. Export isolate-level PCoA coordinates
# ------------------------------------------------------------

write.csv(
  pcoa_df,
  "PCoA_coordinates_population.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Export publication-quality PCoA figure
# ------------------------------------------------------------

ggsave(
  filename = "PCoA_population_structure.tiff",
  plot = pcoa_plot,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


# ------------------------------------------------------------
# 11. Report summary
# ------------------------------------------------------------

cat(
  "\nPCoA analysis completed successfully.\n\n"
)

cat(
  "Number of isolates:",
  nrow(pcoa_df),
  "\n"
)

cat(
  "Number of populations:",
  n_distinct(pcoa_df$District),
  "\n"
)

cat(
  "PCoA1 variance explained:",
  round(pcoa_variance[1], 2),
  "%\n"
)

cat(
  "PCoA2 variance explained:",
  round(pcoa_variance[2], 2),
  "%\n"
)

cat(
  "\nPopulation centroids:\n"
)

print(
  pcoa_centroids
)

print(
  pcoa_plot
)