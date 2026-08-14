# ============================================================
# Discriminant Analysis of Principal Components (DAPC)
#
# Study: Puccinia striiformis f. sp. tritici populations from
#        Khyber Pakhtunkhwa, Pakistan
#
# Analysis:
#   1. Infer genetic clusters using BIC
#   2. Optimize the number of retained principal components
#      using the a-score criterion
#   3. Perform final DAPC
#   4. Visualize genetic clusters
#   5. Examine cluster distribution among districts
#
# Required input objects:
#   gen              : adegenet genind object
#   mlg_geo_correct  : isolate metadata containing Isolate and
#                      District columns
#
# Outputs:
#   DAPC_BIC_cluster_selection.tiff
#   DAPC_population_structure.tiff
#   DAPC_cluster_composition_barplot.tiff
#   DAPC_coordinates.csv
#   DAPC_cluster_by_district.csv
#   DAPC_cluster_composition.csv
#
# ============================================================


# ============================================================
# 1. Load required packages
# ============================================================

library(adegenet)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)


# ============================================================
# 2. Validate required input objects
# ============================================================

required_objects <- c(
  "gen",
  "mlg_geo_correct"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1))
]

if (length(missing_objects) > 0) {
  
  stop(
    "The following required objects are missing: ",
    paste(missing_objects, collapse = ", ")
  )
}


# Check required metadata columns

required_columns <- c(
  "Isolate",
  "District"
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
# PART I: CLUSTER SELECTION USING BIC
# ============================================================


# ------------------------------------------------------------
# 3. Evaluate candidate numbers of genetic clusters
# ------------------------------------------------------------

set.seed(123)

bic_test <- find.clusters(
  gen,
  max.n.clust = 10,
  n.pca = 50,
  choose.n.clust = FALSE
)


# Extract BIC values

bic_df <- data.frame(
  K = seq_along(bic_test$Kstat),
  BIC = as.numeric(bic_test$Kstat)
)


# ------------------------------------------------------------
# 4. Selected number of clusters
# ------------------------------------------------------------

selected_k <- 6


# ------------------------------------------------------------
# 5. Create BIC plot
# ------------------------------------------------------------

bic_plot <- ggplot(
  bic_df,
  aes(
    x = K,
    y = BIC
  )
) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_point(
    size = 3
  ) +
  
  geom_point(
    data = subset(
      bic_df,
      K == selected_k
    ),
    size = 5
  ) +
  
  geom_vline(
    xintercept = selected_k,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  annotate(
    "text",
    x = selected_k,
    y = max(bic_df$BIC) * 0.85,
    label = paste0(
      "K = ",
      selected_k
    ),
    size = 5,
    fontface = "bold",
    hjust = -0.2
  ) +
  
  scale_x_continuous(
    breaks = seq_len(max(bic_df$K))
  ) +
  
  labs(
    x = "Number of clusters (K)",
    y = "Bayesian Information Criterion (BIC)"
  ) +
  
  theme_classic(
    base_size = 15
  ) +
  
  theme(
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    axis.ticks.length = grid::unit(
      0.15,
      "cm"
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 15
    ),
    
    axis.text = element_text(
      size = 12
    )
  )


# Export BIC plot

ggsave(
  filename = "DAPC_BIC_cluster_selection.tiff",
  plot = bic_plot,
  width = 7,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


# ============================================================
# PART II: ASSIGN GENETIC CLUSTERS
# ============================================================


set.seed(123)

dapc_clusters <- find.clusters(
  gen,
  max.n.clust = 10,
  n.pca = 50,
  n.clust = selected_k
)


# ============================================================
# PART III: OPTIMIZE NUMBER OF RETAINED PCs
# ============================================================


# ------------------------------------------------------------
# 6. Initial DAPC
# ------------------------------------------------------------

set.seed(123)

dapc_initial <- dapc(
  gen,
  pop = dapc_clusters$grp,
  n.pca = 50,
  n.da = selected_k - 1
)


# ------------------------------------------------------------
# 7. Optimize retained principal components using a-score
# ------------------------------------------------------------

set.seed(123)

optim_pca <- optim.a.score(
  dapc_initial,
  n.pca = seq(10, 44, by = 2),
  n.da = selected_k - 1
)


# ------------------------------------------------------------
# 8. Selected number of principal components
# ------------------------------------------------------------

selected_n_pca <- 20


# ============================================================
# PART IV: FINAL DAPC
# ============================================================


selected_n_da <- selected_k - 1


set.seed(123)

dapc_final <- dapc(
  gen,
  pop = dapc_clusters$grp,
  n.pca = selected_n_pca,
  n.da = selected_n_da
)


# ============================================================
# PART V: EXTRACT DAPC COORDINATES
# ============================================================


dapc_coords <- as.data.frame(
  dapc_final$ind.coord[, 1:2]
)

colnames(dapc_coords) <- c(
  "LD1",
  "LD2"
)


# Add isolate identifiers

dapc_coords$Isolate <- rownames(
  dapc_coords
)


# Add cluster assignments

dapc_coords$Cluster <- factor(
  dapc_final$grp,
  levels = seq_len(selected_k),
  labels = paste0(
    "Cluster ",
    seq_len(selected_k)
  )
)


# ============================================================
# PART VI: CALCULATE CLUSTER CENTROIDS
# ============================================================


dapc_centroids <- dapc_coords %>%
  group_by(Cluster) %>%
  summarise(
    LD1 = mean(LD1),
    LD2 = mean(LD2),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Cluster_label = paste0(
      Cluster,
      " (n = ",
      n,
      ")"
    )
  )


# ============================================================
# PART VII: DEFINE CLUSTER COLOURS
# ============================================================


dapc_colors <- c(
  "Cluster 1" = "#1B9E77",
  "Cluster 2" = "#D95F02",
  "Cluster 3" = "#7570B3",
  "Cluster 4" = "#E7298A",
  "Cluster 5" = "#66A61E",
  "Cluster 6" = "#E6AB02"
)


# ============================================================
# PART VIII: PREPARE DATA FOR CONFIDENCE ELLIPSES
# ============================================================


ellipse_data <- dapc_coords %>%
  group_by(Cluster) %>%
  filter(
    sd(LD1) > 0,
    sd(LD2) > 0
  ) %>%
  ungroup()


# ============================================================
# PART IX: CREATE FINAL DAPC PLOT
# ============================================================


dapc_plot_final <- ggplot(
  dapc_coords,
  aes(
    x = LD1,
    y = LD2,
    color = Cluster
  )
) +
  
  # Reference axes
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
  
  # Cluster ellipses
  stat_ellipse(
    data = ellipse_data,
    aes(
      fill = Cluster
    ),
    geom = "polygon",
    alpha = 0.15,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  # Individual isolates
  geom_point(
    size = 3,
    alpha = 0.75
  ) +
  
  # Cluster centroids
  geom_point(
    data = dapc_centroids,
    aes(
      x = LD1,
      y = LD2,
      color = Cluster
    ),
    shape = 4,
    size = 6,
    stroke = 1.5
  ) +
  
  # Cluster labels
  geom_text_repel(
    data = dapc_centroids,
    aes(
      x = LD1,
      y = LD2,
      label = Cluster_label,
      color = Cluster
    ),
    size = 4.5,
    fontface = "bold",
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = Inf,
    segment.color = "grey50"
  ) +
  
  scale_color_manual(
    values = dapc_colors,
    drop = FALSE
  ) +
  
  scale_fill_manual(
    values = dapc_colors,
    drop = FALSE
  ) +
  
  labs(
    x = "Discriminant function 1 (LD1)",
    y = "Discriminant function 2 (LD2)",
    color = "Genetic cluster"
  ) +
  
  theme_classic(
    base_size = 15
  ) +
  
  theme(
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    axis.ticks.length = grid::unit(
      0.15,
      "cm"
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 15
    ),
    
    axis.text = element_text(
      size = 12
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 13
    ),
    
    legend.text = element_text(
      size = 11
    ),
    
    legend.position = "right"
  )


# Export DAPC plot

ggsave(
  filename = "DAPC_population_structure.tiff",
  plot = dapc_plot_final,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


# ============================================================
# PART X: CLUSTER DISTRIBUTION AMONG DISTRICTS
# ============================================================


dapc_cluster_geo <- dapc_coords %>%
  left_join(
    mlg_geo_correct %>%
      select(
        Isolate,
        District
      ),
    by = "Isolate"
  )


# Check for unmatched isolates

if (any(is.na(dapc_cluster_geo$District))) {
  
  warning(
    "Some isolates could not be matched to district metadata."
  )
}


# Create cluster-by-district table

cluster_by_district <- table(
  dapc_cluster_geo$Cluster,
  dapc_cluster_geo$District
)


# ============================================================
# PART XI: CALCULATE DISTRICT-LEVEL CLUSTER COMPOSITION
# ============================================================


cluster_comp <- dapc_cluster_geo %>%
  count(
    District,
    Cluster
  ) %>%
  complete(
    District,
    Cluster,
    fill = list(n = 0)
  ) %>%
  group_by(District) %>%
  mutate(
    Percent = n / sum(n) * 100
  ) %>%
  ungroup()


# Define manuscript district order

district_order <- c(
  "Peshawar",
  "Kohat",
  "Charsadda",
  "Mardan",
  "Swabi",
  "Manshera",
  "Buner"
)


cluster_comp$District <- factor(
  cluster_comp$District,
  levels = district_order
)


cluster_comp$Cluster <- factor(
  cluster_comp$Cluster,
  levels = paste0(
    "Cluster ",
    seq_len(selected_k)
  )
)


# District sample sizes

district_n <- cluster_comp %>%
  group_by(District) %>%
  summarise(
    n = sum(n),
    .groups = "drop"
  )


# ============================================================
# PART XII: CREATE CLUSTER COMPOSITION BAR PLOT
# ============================================================


cluster_barplot_final <- ggplot(
  cluster_comp,
  aes(
    x = District,
    y = Percent,
    fill = Cluster
  )
) +
  
  geom_bar(
    stat = "identity",
    width = 0.65,
    color = NA
  ) +
  
  geom_text(
    data = district_n,
    aes(
      x = District,
      y = 103,
      label = paste0(
        "n = ",
        n
      )
    ),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = dapc_colors
  ) +
  
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 110),
    breaks = seq(0, 100, 25)
  ) +
  
  labs(
    x = NULL,
    y = "Genetic cluster composition (%)",
    fill = "DAPC cluster"
  ) +
  
  theme_classic(
    base_size = 15
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      size = 13,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      face = "bold",
      size = 15
    ),
    
    legend.position = "bottom",
    
    legend.title = element_text(
      face = "bold",
      size = 13
    ),
    
    legend.text = element_text(
      size = 12
    ),
    
    legend.key.height = grid::unit(
      0.4,
      "cm"
    ),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    )
  )


# Export cluster composition plot

ggsave(
  filename = "DAPC_cluster_composition_barplot.tiff",
  plot = cluster_barplot_final,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


# ============================================================
# PART XIII: EXPORT NUMERICAL RESULTS
# ============================================================


# Individual DAPC coordinates and cluster assignments

write.csv(
  dapc_coords,
  "DAPC_coordinates.csv",
  row.names = FALSE
)


# Cluster-by-district counts

write.csv(
  as.data.frame.matrix(cluster_by_district),
  "DAPC_cluster_by_district.csv",
  row.names = TRUE
)


# District-level cluster composition

write.csv(
  cluster_comp,
  "DAPC_cluster_composition.csv",
  row.names = FALSE
)


# ============================================================
# PART XIV: REPORT SUMMARY
# ============================================================


cat(
  "\nDAPC analysis completed successfully.\n\n"
)

cat(
  "Selected number of genetic clusters:",
  selected_k,
  "\n"
)

cat(
  "Retained principal components:",
  selected_n_pca,
  "\n"
)

cat(
  "Retained discriminant functions:",
  selected_n_da,
  "\n\n"
)

cat(
  "Cluster sizes:\n"
)

print(
  table(dapc_coords$Cluster)
)

cat(
  "\nCluster distribution among districts:\n"
)

print(
  cluster_by_district
)