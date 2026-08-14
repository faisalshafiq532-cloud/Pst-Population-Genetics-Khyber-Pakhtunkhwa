# ============================================================
# Spatial Visualization of Estimated Gene Flow Among Populations
#
# Study: Puccinia striiformis f. sp. tritici populations from
#        Khyber Pakhtunkhwa, Pakistan
#
# Required input objects:
#   mlg_geo      : isolate-level geographic and population data
#   Nm_final     : pairwise estimated gene flow matrix
#
# Analyses:
#   1. Calculate geographic centroids of district populations
#   2. Convert pairwise Nm estimates into a network edge list
#   3. Select the strongest estimated gene-flow connections
#   4. Generate a spatial gene-flow network over a terrain map
#
# Output:
#   Gene_Flow_Network/
#     Estimated_Gene_Flow_Network.tiff
#     Estimated_Gene_Flow_Edges.csv
#     Population_Coordinates.csv
#
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(dplyr)
library(sf)
library(terra)
library(geodata)
library(ggplot2)
library(ggspatial)
library(ggnewscale)


# ------------------------------------------------------------
# 2. Define output directory
# ------------------------------------------------------------

output_dir <- "Gene_Flow_Network"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# ============================================================
# PART I: VALIDATE INPUT DATA
# ============================================================


# ------------------------------------------------------------
# 3. Check required objects
# ------------------------------------------------------------

required_objects <- c(
  "mlg_geo",
  "Nm_final"
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
# 4. Verify required geographic columns
# ------------------------------------------------------------

required_columns <- c(
  "District",
  "Longitude (E)",
  "Latitude (N)"
)

missing_columns <- required_columns[
  !required_columns %in% names(mlg_geo)
]

if (length(missing_columns) > 0) {
  
  stop(
    "The following columns are missing from mlg_geo: ",
    paste(missing_columns, collapse = ", ")
  )
}


# ============================================================
# PART II: POPULATION GEOGRAPHIC COORDINATES
# ============================================================


# ------------------------------------------------------------
# 5. Calculate geographic centroids of sampled populations
# ------------------------------------------------------------

population_coords <- mlg_geo %>%
  filter(
    !is.na(District)
  ) %>%
  group_by(District) %>%
  summarise(
    Longitude = mean(
      `Longitude (E)`,
      na.rm = TRUE
    ),
    
    Latitude = mean(
      `Latitude (N)`,
      na.rm = TRUE
    ),
    
    n = n(),
    
    .groups = "drop"
  ) %>%
  rename(
    Population = District
  ) %>%
  mutate(
    Population = as.character(Population)
  )


# ============================================================
# PART III: PREPARE ESTIMATED GENE FLOW NETWORK
# ============================================================


# ------------------------------------------------------------
# 6. Convert Nm matrix to edge list
# ------------------------------------------------------------

Nm_edges <- as.data.frame(
  as.table(Nm_final)
)

colnames(Nm_edges) <- c(
  "From",
  "To",
  "Nm"
)


# ------------------------------------------------------------
# 7. Remove self-comparisons
# ------------------------------------------------------------

Nm_edges <- Nm_edges %>%
  filter(
    From != To
  )


# ------------------------------------------------------------
# 8. Remove duplicate undirected connections
# ------------------------------------------------------------

Nm_edges <- Nm_edges %>%
  rowwise() %>%
  mutate(
    Pair = paste(
      sort(c(From, To)),
      collapse = "_"
    )
  ) %>%
  ungroup() %>%
  distinct(
    Pair,
    .keep_all = TRUE
  ) %>%
  select(
    -Pair
  )


# ------------------------------------------------------------
# 9. Convert Nm values to numeric
#
# Values reported as >10 are assigned a plotting value of 11
# solely for visualization and ranking.
# ------------------------------------------------------------

Nm_edges <- Nm_edges %>%
  mutate(
    Nm = as.character(Nm),
    
    Nm_numeric = if_else(
      Nm == ">10",
      11,
      as.numeric(Nm)
    )
  )


# ------------------------------------------------------------
# 10. Select strongest estimated gene-flow connections
# ------------------------------------------------------------

Nm_top_connections <- Nm_edges %>%
  arrange(
    desc(Nm_numeric)
  ) %>%
  slice_head(
    n = min(12, n())
  )


# ============================================================
# PART IV: ADD GEOGRAPHIC COORDINATES TO NETWORK EDGES
# ============================================================


# ------------------------------------------------------------
# 11. Add coordinates of source and destination populations
# ------------------------------------------------------------

Nm_edges_spatial <- Nm_top_connections %>%
  
  left_join(
    population_coords,
    by = c("From" = "Population")
  ) %>%
  
  rename(
    From_lon = Longitude,
    From_lat = Latitude,
    From_n = n
  ) %>%
  
  left_join(
    population_coords,
    by = c("To" = "Population")
  ) %>%
  
  rename(
    To_lon = Longitude,
    To_lat = Latitude,
    To_n = n
  )


# ------------------------------------------------------------
# 12. Check for failed population-coordinate matches
# ------------------------------------------------------------

if (
  any(
    is.na(Nm_edges_spatial$From_lon)
  ) ||
  any(
    is.na(Nm_edges_spatial$To_lon)
  )
) {
  
  stop(
    "Some populations in Nm_final could not be matched ",
    "to geographic coordinates."
  )
}


# ============================================================
# PART V: CREATE SPATIAL LINE GEOMETRIES
# ============================================================


# ------------------------------------------------------------
# 13. Convert population connections to spatial lines
# ------------------------------------------------------------

Nm_lines <- lapply(
  seq_len(nrow(Nm_edges_spatial)),
  
  function(i) {
    
    st_linestring(
      matrix(
        c(
          Nm_edges_spatial$From_lon[i],
          Nm_edges_spatial$From_lat[i],
          Nm_edges_spatial$To_lon[i],
          Nm_edges_spatial$To_lat[i]
        ),
        
        ncol = 2,
        byrow = TRUE
      )
    )
  }
)


# ------------------------------------------------------------
# 14. Create spatial network object
# ------------------------------------------------------------

Nm_lines_sf <- st_sf(
  Nm_edges_spatial,
  
  geometry = st_sfc(
    Nm_lines,
    crs = 4326
  )
)


# ------------------------------------------------------------
# 15. Convert population coordinates to spatial points
# ------------------------------------------------------------

population_nodes_sf <- st_as_sf(
  population_coords,
  
  coords = c(
    "Longitude",
    "Latitude"
  ),
  
  crs = 4326
)


# ============================================================
# PART VI: PREPARE KPK BASE MAP
# ============================================================


# ------------------------------------------------------------
# 16. Download Pakistan administrative boundaries
# ------------------------------------------------------------

pak_admin3 <- geodata::gadm(
  country = "PAK",
  level = 3,
  path = output_dir
)


# ------------------------------------------------------------
# 17. Extract Khyber Pakhtunkhwa boundaries
# ------------------------------------------------------------

kpk_districts <- pak_admin3[
  pak_admin3$NAME_1 == "Khyber-Pakhtunkhwa",
]


# ------------------------------------------------------------
# 18. Convert boundaries to sf format
# ------------------------------------------------------------

kpk_districts_sf <- st_as_sf(
  kpk_districts
)


# ============================================================
# PART VII: CREATE TERRAIN BACKGROUND
# ============================================================


# ------------------------------------------------------------
# 19. Download elevation data
# ------------------------------------------------------------

pak_dem <- geodata::elevation_30s(
  country = "PAK",
  path = output_dir
)


# ------------------------------------------------------------
# 20. Crop elevation data to KPK region
# ------------------------------------------------------------

kpk_extent <- terra::ext(
  70.5,
  74.2,
  33.0,
  36.0
)

kpk_dem <- terra::crop(
  pak_dem,
  kpk_extent
)


# ------------------------------------------------------------
# 21. Calculate hillshade
# ------------------------------------------------------------

slope <- terra::terrain(
  kpk_dem,
  v = "slope",
  unit = "radians"
)

aspect <- terra::terrain(
  kpk_dem,
  v = "aspect",
  unit = "radians"
)

kpk_hillshade <- terra::shade(
  slope,
  aspect,
  angle = 45,
  direction = 225
)


# ------------------------------------------------------------
# 22. Convert hillshade to data frame
# ------------------------------------------------------------

hillshade_df <- as.data.frame(
  kpk_hillshade,
  xy = TRUE
)

colnames(hillshade_df)[3] <- "hillshade"


# ============================================================
# PART VIII: DEFINE POPULATION COLOURS
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


population_nodes_sf$color <- unname(
  district_colors[
    population_nodes_sf$Population
  ]
)


# ============================================================
# PART IX: CREATE ESTIMATED GENE FLOW NETWORK MAP
# ============================================================


gene_flow_map <- ggplot() +
  
  # Terrain background
  geom_raster(
    data = hillshade_df,
    
    aes(
      x = x,
      y = y,
      fill = hillshade
    ),
    
    alpha = 0.55
  ) +
  
  scale_fill_gradient(
    low = "grey95",
    high = "grey30",
    guide = "none"
  ) +
  
  new_scale_fill() +
  
  # District boundaries
  geom_sf(
    data = kpk_districts_sf,
    fill = NA,
    color = "grey60",
    linewidth = 0.25
  ) +
  
  # Estimated gene-flow connections
  geom_sf(
    data = Nm_lines_sf,
    
    aes(
      linewidth = log10(Nm_numeric + 1)
    ),
    
    color = "red",
    alpha = 0.65
  ) +
  
  scale_linewidth_continuous(
    name = "Estimated gene flow\n(log10 Nm + 1)",
    range = c(0.4, 2)
  ) +
  
  # Population nodes
  geom_sf(
    data = population_nodes_sf,
    
    aes(
      size = n,
      fill = color
    ),
    
    shape = 21,
    color = "black",
    stroke = 1
  ) +
  
  scale_fill_identity() +
  
  scale_size_continuous(
    name = "Number of isolates",
    range = c(3, 8)
  ) +
  
  # Population labels
  geom_sf_text(
    data = population_nodes_sf,
    
    aes(
      label = Population
    ),
    
    nudge_y = 0.05,
    size = 4,
    fontface = "bold"
  ) +
  
  coord_sf(
    xlim = c(71, 74),
    ylim = c(33.4, 34.7),
    expand = FALSE
  ) +
  
  annotation_scale(
    location = "bl",
    width_hint = 0.25
  ) +
  
  annotation_north_arrow(
    location = "tl",
    which_north = "true",
    height = unit(0.8, "cm"),
    width = unit(0.8, "cm")
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    axis.title = element_blank(),
    
    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = 0.8
    ),
    
    axis.ticks = element_line(
      color = "black"
    ),
    
    legend.position = "right"
  )


# ============================================================
# PART X: EXPORT RESULTS
# ============================================================


# ------------------------------------------------------------
# 23. Export selected estimated gene-flow connections
# ------------------------------------------------------------

write.csv(
  Nm_edges_spatial,
  
  file = file.path(
    output_dir,
    "Estimated_Gene_Flow_Edges.csv"
  ),
  
  row.names = FALSE
)


# ------------------------------------------------------------
# 24. Export population geographic coordinates
# ------------------------------------------------------------

write.csv(
  population_coords,
  
  file = file.path(
    output_dir,
    "Population_Coordinates.csv"
  ),
  
  row.names = FALSE
)


# ------------------------------------------------------------
# 25. Export publication-quality figure
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "Estimated_Gene_Flow_Network.tiff"
  ),
  
  plot = gene_flow_map,
  
  width = 230,
  height = 160,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)


# ------------------------------------------------------------
# 26. Report summary
# ------------------------------------------------------------

cat(
  "\nSpatial visualization of estimated gene flow complete.\n\n"
)

cat(
  "Number of populations:",
  nrow(population_coords),
  "\n"
)

cat(
  "Total pairwise connections:",
  nrow(Nm_edges),
  "\n"
)

cat(
  "Connections displayed:",
  nrow(Nm_edges_spatial),
  "\n"
)

cat(
  "\nOutput directory:",
  output_dir,
  "\n"
)

print(
  Nm_edges_spatial
)

print(
  gene_flow_map
)