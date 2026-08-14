# ============================================================
# Atmospheric Connectivity Analysis
#
# Study:
# Puccinia striiformis f. sp. tritici populations from
# Khyber Pakhtunkhwa, Pakistan
#
# Purpose:
# Evaluate potential wind-mediated connectivity among seven
# sampled populations and test its association with genetic
# differentiation.
#
# Data availability:
# The large ERA5 hourly wind dataset and geographic boundary
# files used in this analysis are not deposited in this GitHub
# repository.
#
# ERA5 data:
# Hourly 10-m zonal (u) and meridional (v) wind components
# for February-May during 2000-2025 were obtained separately
# and extracted at population centroid coordinates using
# bilinear interpolation.
#
# To reproduce this analysis, users must obtain the relevant
# ERA5 data and prepare the following objects before running
# the script:
#
#   wind_u_pop         Hourly zonal wind component matrix
#   wind_v_pop         Hourly meridional wind component matrix
#   population_coords  Population centroid coordinates with:
#                      District, Longitude, and Latitude
#   pairwise_fst       Pairwise Weir-Cockerham FST matrix
#   kpk_districts_sf   Khyber Pakhtunkhwa district boundaries
#                      as an sf object
#
# Rows of wind_u_pop and wind_v_pop represent identical hourly
# observations, and columns represent the same populations in
# the same order.
#
# This script calculates:
#   1. Wind speed and direction
#   2. Population-to-population geographic bearings
#   3. Hourly directional alignment
#   4. Wind-speed-weighted atmospheric transport
#   5. Atmospheric Connectivity Index (ACI)
#   6. Symmetric atmospheric connectivity
#   7. Dominant atmospheric pathways
#   8. Atmospheric connectivity heatmap
#   9. Wind rose summaries
#  10. Mantel association between atmospheric isolation and
#      genetic differentiation
#
# The final publication figures are provided in the manuscript
# and are not exported by this script.
# ============================================================


# ============================================================
# 1. Load required packages
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(geosphere)
library(vegan)
library(grid)


# ============================================================
# 2. Define sampled populations
# ============================================================

population_order <- c(
  "Peshawar",
  "Kohat",
  "Charsadda",
  "Mardan",
  "Swabi",
  "Mansehra",
  "Buner"
)


# ============================================================
# 3. Validate required objects
# ============================================================

required_objects <- c(
  "wind_u_pop",
  "wind_v_pop",
  "population_coords",
  "pairwise_fst",
  "kpk_districts_sf"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    inherits = FALSE
  )
]

if (length(missing_objects) > 0) {
  
  stop(
    "The following required objects are missing:\n",
    paste(
      missing_objects,
      collapse = "\n"
    )
  )
}


# ============================================================
# 4. Validate population coordinate data
# ============================================================

required_coordinate_columns <- c(
  "District",
  "Longitude",
  "Latitude"
)

missing_coordinate_columns <- required_coordinate_columns[
  !required_coordinate_columns %in%
    names(population_coords)
]

if (length(missing_coordinate_columns) > 0) {
  
  stop(
    "population_coords is missing the following columns: ",
    paste(
      missing_coordinate_columns,
      collapse = ", "
    )
  )
}


population_coords <- population_coords %>%
  filter(
    District %in% population_order
  ) %>%
  mutate(
    District = factor(
      District,
      levels = population_order
    )
  ) %>%
  arrange(District)


if (nrow(population_coords) != length(population_order)) {
  
  stop(
    "population_coords must contain one row for each ",
    "of the seven sampled populations."
  )
}


# ============================================================
# 5. Validate and standardize population names
# ============================================================

colnames(wind_u_pop) <- trimws(
  colnames(wind_u_pop)
)

colnames(wind_v_pop) <- trimws(
  colnames(wind_v_pop)
)

rownames(pairwise_fst) <- trimws(
  rownames(pairwise_fst)
)

colnames(pairwise_fst) <- trimws(
  colnames(pairwise_fst)
)


# Standardize the spelling of Mansehra where necessary

colnames(wind_u_pop)[
  colnames(wind_u_pop) == "Manshera"
] <- "Mansehra"

colnames(wind_v_pop)[
  colnames(wind_v_pop) == "Manshera"
] <- "Mansehra"

rownames(pairwise_fst)[
  rownames(pairwise_fst) == "Manshera"
] <- "Mansehra"

colnames(pairwise_fst)[
  colnames(pairwise_fst) == "Manshera"
] <- "Mansehra"


# ============================================================
# 6. Validate wind matrices
# ============================================================

if (!identical(
  dim(wind_u_pop),
  dim(wind_v_pop)
)) {
  
  stop(
    "wind_u_pop and wind_v_pop must have identical dimensions."
  )
}


if (!identical(
  colnames(wind_u_pop),
  colnames(wind_v_pop)
)) {
  
  stop(
    "wind_u_pop and wind_v_pop must have identical ",
    "population columns in the same order."
  )
}


if (!all(
  population_order %in%
  colnames(wind_u_pop)
)) {
  
  stop(
    "One or more sampled populations are missing ",
    "from the wind data."
  )
}


wind_u_pop <- wind_u_pop[
  ,
  population_order,
  drop = FALSE
]

wind_v_pop <- wind_v_pop[
  ,
  population_order,
  drop = FALSE
]


# ============================================================
# 7. Calculate wind speed
# ============================================================

wind_speed_pop <- sqrt(
  wind_u_pop^2 +
    wind_v_pop^2
)


# ============================================================
# 8. Calculate wind direction
#
# Meteorological convention:
# Direction represents the direction toward which the wind is
# moving, allowing direct comparison with source-to-destination
# population bearings.
# ============================================================

wind_dir_pop <- (
  atan2(
    wind_u_pop,
    wind_v_pop
  ) *
    180 / pi +
    360
) %% 360


# ============================================================
# 9. Calculate pairwise geographic distances
# ============================================================

population_lonlat <- population_coords %>%
  select(
    Longitude,
    Latitude
  ) %>%
  as.matrix()


dist_matrix <- geosphere::distm(
  population_lonlat,
  fun = geosphere::distHaversine
) / 1000


rownames(dist_matrix) <- population_order
colnames(dist_matrix) <- population_order


diag(dist_matrix) <- NA


# ============================================================
# 10. Calculate source-to-destination bearings
# ============================================================

n_pop <- length(population_order)

bearing_matrix <- matrix(
  NA_real_,
  nrow = n_pop,
  ncol = n_pop,
  dimnames = list(
    population_order,
    population_order
  )
)


for (i in seq_len(n_pop)) {
  
  for (j in seq_len(n_pop)) {
    
    if (i != j) {
      
      bearing_matrix[i, j] <- geosphere::bearing(
        p1 = population_lonlat[i, ],
        p2 = population_lonlat[j, ]
      )
      
      bearing_matrix[i, j] <- (
        bearing_matrix[i, j] +
          360
      ) %% 360
    }
  }
}


# ============================================================
# 11. Define directional alignment function
#
# Values range from:
#   1 = complete directional agreement
#   0 = opposite direction
# ============================================================

wind_alignment <- function(
    wind_direction,
    target_direction
) {
  
  angle_difference <- abs(
    (
      (
        wind_direction -
          target_direction +
          180
      ) %% 360
    ) - 180
  )
  
  (
    1 +
      cos(
        angle_difference *
          pi / 180
      )
  ) / 2
}


# ============================================================
# 12. Calculate hourly directional alignment
# ============================================================

alignment_hourly <- vector(
  "list",
  nrow(wind_dir_pop)
)


for (t in seq_len(nrow(wind_dir_pop))) {
  
  alignment_matrix <- matrix(
    NA_real_,
    nrow = n_pop,
    ncol = n_pop,
    dimnames = list(
      population_order,
      population_order
    )
  )
  
  for (i in seq_len(n_pop)) {
    
    for (j in seq_len(n_pop)) {
      
      if (i != j) {
        
        alignment_matrix[i, j] <- wind_alignment(
          wind_direction = wind_dir_pop[t, i],
          target_direction = bearing_matrix[i, j]
        )
      }
    }
  }
  
  alignment_hourly[[t]] <- alignment_matrix
}


# ============================================================
# 13. Calculate hourly wind transport
#
# Directional alignment is weighted by wind speed at the
# source population.
# ============================================================

transport_hourly <- vector(
  "list",
  nrow(wind_speed_pop)
)


for (t in seq_len(nrow(wind_speed_pop))) {
  
  transport_matrix <- alignment_hourly[[t]]
  
  for (i in seq_len(n_pop)) {
    
    transport_matrix[i, ] <-
      transport_matrix[i, ] *
      wind_speed_pop[t, i]
  }
  
  diag(transport_matrix) <- NA
  
  transport_hourly[[t]] <- transport_matrix
}


# ============================================================
# 14. Calculate mean atmospheric transport
# ============================================================

transport_mean <- matrix(
  NA_real_,
  nrow = n_pop,
  ncol = n_pop,
  dimnames = list(
    population_order,
    population_order
  )
)


for (i in seq_len(n_pop)) {
  
  for (j in seq_len(n_pop)) {
    
    if (i != j) {
      
      transport_mean[i, j] <- mean(
        vapply(
          transport_hourly,
          function(x) x[i, j],
          numeric(1)
        ),
        na.rm = TRUE
      )
    }
  }
}


diag(transport_mean) <- 0


# ============================================================
# 15. Calculate Atmospheric Connectivity Index (ACI)
#
# Higher values indicate greater directional and wind-speed-
# weighted connectivity relative to geographic separation.
# ============================================================

ACI <- transport_mean /
  dist_matrix


diag(ACI) <- 0

ACI[
  is.infinite(ACI)
] <- 0

ACI[
  is.na(ACI)
] <- 0


# ============================================================
# 16. Calculate symmetric atmospheric connectivity
#
# Reciprocal directional connections are averaged for analyses
# requiring an undirected pairwise matrix.
# ============================================================

ACI_sym <- (
  ACI +
    t(ACI)
) / 2


diag(ACI_sym) <- 0

rownames(ACI_sym) <- population_order
colnames(ACI_sym) <- population_order


# ============================================================
# 17. Create directed atmospheric connectivity edge table
# ============================================================

ACI_edges <- expand.grid(
  Source = rownames(ACI),
  Destination = colnames(ACI),
  stringsAsFactors = FALSE
) %>%
  mutate(
    ACI = as.vector(ACI)
  ) %>%
  filter(
    Source != Destination
  ) %>%
  left_join(
    population_coords %>%
      mutate(
        District = as.character(District)
      ) %>%
      select(
        District,
        Longitude,
        Latitude
      ),
    by = c(
      "Source" = "District"
    )
  ) %>%
  rename(
    lon_source = Longitude,
    lat_source = Latitude
  ) %>%
  left_join(
    population_coords %>%
      mutate(
        District = as.character(District)
      ) %>%
      select(
        District,
        Longitude,
        Latitude
      ),
    by = c(
      "Destination" = "District"
    )
  ) %>%
  rename(
    lon_dest = Longitude,
    lat_dest = Latitude
  )


# ============================================================
# 18. Select dominant atmospheric pathways
# ============================================================

n_display_pathways <- 15


ACI_plot_edges <- ACI_edges %>%
  arrange(
    desc(ACI)
  ) %>%
  slice_head(
    n = n_display_pathways
  ) %>%
  mutate(
    ACI_scaled = ACI / max(ACI)
  )


# ============================================================
# 19. Prepare Khyber Pakhtunkhwa map
# ============================================================

sampled_districts <- population_order


sampled_kpk_sf <- kpk_districts_sf %>%
  filter(
    NAME_3 %in% sampled_districts
  )


sampled_kpk_dissolved <- sampled_kpk_sf %>%
  group_by(NAME_3) %>%
  summarise(
    geometry = sf::st_union(geometry),
    .groups = "drop"
  )


# ============================================================
# 20. Create atmospheric connectivity network
# ============================================================

kpk_base_map <- ggplot() +
  
  geom_sf(
    data = kpk_districts_sf,
    fill = "grey95",
    color = "grey75",
    linewidth = 0.25
  ) +
  
  geom_sf(
    data = sampled_kpk_dissolved,
    fill = "grey80",
    color = "black",
    linewidth = 0.5
  ) +
  
  geom_point(
    data = population_coords,
    aes(
      x = Longitude,
      y = Latitude
    ),
    size = 3
  ) +
  
  coord_sf(
    xlim = c(71.0, 74.3),
    ylim = c(33.0, 35.3),
    expand = FALSE
  ) +
  
  labs(
    x = "Longitude (°E)",
    y = "Latitude (°N)"
  ) +
  
  theme_classic()


kpk_aci_map <- kpk_base_map +
  
  geom_curve(
    data = ACI_plot_edges,
    aes(
      x = lon_source,
      y = lat_source,
      xend = lon_dest,
      yend = lat_dest,
      linewidth = ACI_scaled
    ),
    curvature = 0.15,
    arrow = arrow(
      length = unit(
        0.15,
        "cm"
      ),
      type = "closed"
    ),
    alpha = 0.7,
    color = "black"
  ) +
  
  scale_linewidth_continuous(
    range = c(
      0.3,
      1.5
    ),
    guide = "none"
  )


print(kpk_aci_map)


# ============================================================
# 21. Prepare ACI heatmap data
# ============================================================

ACI_heatmap_df <- as.data.frame(
  ACI_sym
) %>%
  mutate(
    Source = rownames(.)
  ) %>%
  pivot_longer(
    cols = -Source,
    names_to = "Destination",
    values_to = "ACI"
  ) %>%
  mutate(
    Source = factor(
      Source,
      levels = population_order
    ),
    Destination = factor(
      Destination,
      levels = population_order
    ),
    ACI_plot = ifelse(
      Source == Destination,
      NA,
      ACI
    )
  )


# ============================================================
# 22. Plot atmospheric connectivity heatmap
# ============================================================

ACI_heatmap <- ggplot(
  ACI_heatmap_df,
  aes(
    x = Destination,
    y = Source,
    fill = ACI_plot
  )
) +
  
  geom_tile(
    color = "white"
  ) +
  
  scale_fill_viridis_c(
    option = "magma",
    na.value = "white",
    name = "ACI"
  ) +
  
  coord_equal() +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_classic() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


print(ACI_heatmap)


# ============================================================
# 23. Prepare wind direction classes
# ============================================================

wind_direction_classes <- c(
  "N",
  "NE",
  "E",
  "SE",
  "S",
  "SW",
  "W",
  "NW"
)


classify_wind_direction <- function(direction) {
  
  cut(
    direction,
    breaks = c(
      -22.5,
      22.5,
      67.5,
      112.5,
      157.5,
      202.5,
      247.5,
      292.5,
      337.5,
      382.5
    ),
    labels = c(
      "N",
      "NE",
      "E",
      "SE",
      "S",
      "SW",
      "W",
      "NW",
      "N"
    ),
    include.lowest = TRUE,
    right = FALSE
  )
}


# ============================================================
# 24. Create wind summary data
# ============================================================

wind_summary <- lapply(
  seq_len(n_pop),
  function(i) {
    
    data.frame(
      Population = population_order[i],
      Direction = wind_dir_pop[, i],
      Wind_speed = wind_speed_pop[, i]
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    Direction_class = classify_wind_direction(
      Direction
    )
  )


wind_rose_summary <- wind_summary %>%
  filter(
    !is.na(Direction_class)
  ) %>%
  group_by(
    Population,
    Direction_class
  ) %>%
  summarise(
    frequency = 100 * n() / sum(n()),
    mean_speed = mean(
      Wind_speed,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Population = factor(
      Population,
      levels = population_order
    ),
    Direction_class = factor(
      Direction_class,
      levels = wind_direction_classes
    )
  )


# ============================================================
# 25. Plot wind rose summaries
# ============================================================

wind_rose_plot <- ggplot(
  wind_rose_summary,
  aes(
    x = Direction_class,
    y = frequency,
    fill = mean_speed
  )
) +
  
  geom_col(
    width = 1,
    color = "white"
  ) +
  
  coord_polar(
    start = -pi / 8
  ) +
  
  facet_wrap(
    ~Population,
    ncol = 4
  ) +
  
  scale_fill_viridis_c(
    name = "Mean wind\nspeed (m/s)"
  ) +
  
  labs(
    x = NULL,
    y = "Frequency (%)"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text.x = element_text(
      size = 8
    ),
    strip.text = element_text(
      face = "bold"
    )
  )


print(wind_rose_plot)


# ============================================================
# 26. Prepare pairwise FST matrix
# ============================================================

if (!all(
  population_order %in%
  rownames(pairwise_fst)
)) {
  
  stop(
    "One or more sampled populations are missing ",
    "from pairwise_fst."
  )
}


if (!all(
  population_order %in%
  colnames(pairwise_fst)
)) {
  
  stop(
    "One or more sampled populations are missing ",
    "from pairwise_fst."
  )
}


pairwise_fst_ordered <- pairwise_fst[
  population_order,
  population_order
]


# Negative Weir-Cockerham FST estimates are treated as zero
# for the distance-based visualization and Mantel analysis.

fst_matrix_mantel <- pairwise_fst_ordered

fst_matrix_mantel[
  fst_matrix_mantel < 0
] <- 0

diag(fst_matrix_mantel) <- 0


# ============================================================
# 27. Convert atmospheric connectivity to isolation distance
# ============================================================

max_aci <- max(
  ACI_sym,
  na.rm = TRUE
)


if (max_aci <= 0) {
  
  stop(
    "Atmospheric connectivity values are not suitable ",
    "for normalization."
  )
}


ACI_norm <- ACI_sym / max_aci


atmospheric_distance <- 1 -
  ACI_norm


diag(atmospheric_distance) <- 0


# ============================================================
# 28. Perform Mantel test
# ============================================================

fst_dist <- as.dist(
  fst_matrix_mantel
)

atmos_dist <- as.dist(
  atmospheric_distance
)


set.seed(20260814)


mantel_result <- vegan::mantel(
  fst_dist,
  atmos_dist,
  method = "pearson",
  permutations = 9999
)


print(mantel_result)


# ============================================================
# 29. Prepare pairwise data for visualization
# ============================================================

pair_indices <- which(
  upper.tri(fst_matrix_mantel),
  arr.ind = TRUE
)


mantel_df <- data.frame(
  Population1 = rownames(
    fst_matrix_mantel
  )[pair_indices[, 1]],
  
  Population2 = colnames(
    fst_matrix_mantel
  )[pair_indices[, 2]],
  
  FST = fst_matrix_mantel[
    upper.tri(fst_matrix_mantel)
  ],
  
  Atmospheric_distance = atmospheric_distance[
    upper.tri(atmospheric_distance)
  ]
)


# ============================================================
# 30. Plot genetic differentiation against atmospheric isolation
# ============================================================

fst_wind_plot <- ggplot(
  mantel_df,
  aes(
    x = Atmospheric_distance,
    y = FST
  )
) +
  
  geom_point(
    size = 3
  ) +
  
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  
  labs(
    x = "Atmospheric isolation distance (1 - normalized ACI)",
    y = expression(
      "Genetic differentiation (" *
        F[ST] *
        ")"
    )
  ) +
  
  theme_classic() +
  
  theme(
    axis.title = element_text(
      size = 12
    ),
    axis.text = element_text(
      size = 10
    )
  )


print(fst_wind_plot)


# ============================================================
# 31. Report completion
# ============================================================

cat(
  "\nAtmospheric connectivity analysis completed successfully.\n"
)

cat(
  "Populations:",
  n_pop,
  "\n"
)

cat(
  "Wind observations:",
  nrow(wind_u_pop),
  "\n"
)

cat(
  "Dominant atmospheric pathways displayed:",
  nrow(ACI_plot_edges),
  "\n"
)

cat(
  "Mantel statistic:",
  round(
    mantel_result$statistic,
    4
  ),
  "\n"
)

cat(
  "Mantel P-value:",
  mantel_result$signif,
  "\n"
)

cat(
  "Note: ERA5 wind data and geographic boundary files are ",
  "not deposited in this repository. See the script header ",
  "and manuscript Methods section for data preparation details.\n"
)

cat(
  "The final publication figures are provided in the manuscript.\n"
)