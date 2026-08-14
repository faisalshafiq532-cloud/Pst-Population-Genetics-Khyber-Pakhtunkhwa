# ============================================================
# Spatial Distribution of Multilocus Genotypes (MLGs)
#
# Study: Puccinia striiformis f. sp. tritici populations from
#        Khyber Pakhtunkhwa, Pakistan
#
# Inputs:
#   SSR_PK-data.xlsx
#   Locations_Sampling.xlsx
#
# Analyses:
#   1. Create multilocus genotype assignments
#   2. Merge MLG assignments with sampling coordinates
#   3. Identify dominant MLGs (frequency >= 5 isolates)
#   4. Map the geographic distribution of dominant MLGs
#   5. Export isolate-, district-, and MLG-level summaries
#
# Outputs:
#   Spatial_MLG_Distribution/
#     MLG_geographic_distribution_KPK.tiff
#     Supplementary_Table_MLG_distribution_by_district.csv
#     Supplementary_Table_Isolate_MLG_assignment.csv
#     Supplementary_Table_MLG_frequency_summary.csv
#
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(readxl)
library(adegenet)
library(poppr)
library(dplyr)
library(sf)
library(terra)
library(geodata)
library(ggplot2)


# ------------------------------------------------------------
# 2. Define input and output files
# ------------------------------------------------------------

ssr_file <- "SSR_PK-data.xlsx"
coordinate_file <- "Locations_Sampling.xlsx"

output_dir <- "Spatial_MLG_Distribution"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# ============================================================
# PART I: IMPORT AND PREPARE SSR DATA
# ============================================================


# ------------------------------------------------------------
# 3. Import SSR dataset
# ------------------------------------------------------------

ssr_data <- read_excel(
  ssr_file,
  sheet = "PK"
)

# The 17 SSR loci are located in columns 4-20.
locus_cols <- names(ssr_data)[4:20]

ssr_genotypes <- ssr_data[
  locus_cols
]


# ------------------------------------------------------------
# 4. Convert genotype coding to allele-pair format
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
  pop = ssr_data$Regions,
  type = "codom"
)


# ============================================================
# PART II: ASSIGN MULTILOCUS GENOTYPES
# ============================================================


# ------------------------------------------------------------
# 6. Extract isolate-level MLG assignments
# ------------------------------------------------------------

mlg_assign <- poppr::mlg.id(
  gen
)

mlg_df <- data.frame(
  Isolate = unlist(mlg_assign),
  MLG = paste0(
    "MLG_",
    rep(
      names(mlg_assign),
      lengths(mlg_assign)
    )
  ),
  row.names = NULL
)


# ------------------------------------------------------------
# 7. Verify MLG assignment
# ------------------------------------------------------------

stopifnot(
  nrow(mlg_df) == nInd(gen)
)

stopifnot(
  !anyDuplicated(mlg_df$Isolate)
)


# ============================================================
# PART III: MERGE MLG AND SPATIAL INFORMATION
# ============================================================


# ------------------------------------------------------------
# 8. Import sampling coordinates
# ------------------------------------------------------------

coords <- read_excel(
  coordinate_file
)


# ------------------------------------------------------------
# 9. Verify coordinate records
# ------------------------------------------------------------

if (nrow(coords) != nInd(gen)) {
  
  stop(
    "The number of coordinate records does not match ",
    "the number of isolates."
  )
}


# ------------------------------------------------------------
# 10. Assign isolate identifiers
#
# This assumes coordinate records are in the same order as
# isolates in the SSR dataset.
# ------------------------------------------------------------

coords$Isolate <- indNames(
  gen
)


# ------------------------------------------------------------
# 11. Merge MLG assignments and coordinates
# ------------------------------------------------------------

mlg_geo <- coords %>%
  left_join(
    mlg_df,
    by = "Isolate"
  )


# ------------------------------------------------------------
# 12. Validate MLG assignments
# ------------------------------------------------------------

if (any(is.na(mlg_geo$MLG))) {
  
  stop(
    "Some isolates could not be assigned an MLG."
  )
}


# ============================================================
# PART IV: IDENTIFY DOMINANT MLGS
# ============================================================


# ------------------------------------------------------------
# 13. Calculate MLG frequencies
# ------------------------------------------------------------

mlg_frequency <- mlg_geo %>%
  count(
    MLG,
    name = "Number_of_isolates"
  ) %>%
  arrange(
    desc(Number_of_isolates)
  ) %>%
  mutate(
    Frequency_percent = round(
      Number_of_isolates /
        sum(Number_of_isolates) * 100,
      2
    )
  )


# ------------------------------------------------------------
# 14. Identify dominant MLGs
#
# Dominant MLGs are defined as those detected in at least
# five isolates.
# ------------------------------------------------------------

dominant_MLGs <- mlg_frequency %>%
  filter(
    Number_of_isolates >= 5
  ) %>%
  pull(
    MLG
  )


# ------------------------------------------------------------
# 15. Group non-dominant MLGs as "Other"
# ------------------------------------------------------------

mlg_geo <- mlg_geo %>%
  mutate(
    MLG_group = if_else(
      MLG %in% dominant_MLGs,
      MLG,
      "Other"
    )
  )


# ============================================================
# PART V: AGGREGATE MLGS BY SAMPLING LOCATION
# ============================================================


# ------------------------------------------------------------
# 16. Round coordinates
#
# Coordinates are rounded to three decimal places to combine
# isolates collected from the same or nearly identical location.
# ------------------------------------------------------------

mlg_geo <- mlg_geo %>%
  mutate(
    Longitude_round = round(
      `Longitude (E)`,
      3
    ),
    
    Latitude_round = round(
      `Latitude (N)`,
      3
    )
  )


# ------------------------------------------------------------
# 17. Aggregate MLG groups by sampling location
# ------------------------------------------------------------

mlg_map_final <- mlg_geo %>%
  group_by(
    Longitude_round,
    Latitude_round,
    MLG_group
  ) %>%
  summarise(
    N_isolates = n(),
    .groups = "drop"
  )


# ============================================================
# PART VI: OBTAIN KHYBER PAKHTUNKHWA DISTRICT BOUNDARIES
# ============================================================


# ------------------------------------------------------------
# 18. Download administrative boundaries
#
# Pakistan level-3 administrative boundaries are used for
# district-level mapping.
# ------------------------------------------------------------

pak_admin3 <- geodata::gadm(
  country = "PAK",
  level = 3,
  path = output_dir
)


# ------------------------------------------------------------
# 19. Extract Khyber Pakhtunkhwa districts
# ------------------------------------------------------------

kpk_districts <- pak_admin3[
  pak_admin3$NAME_1 == "Khyber-Pakhtunkhwa",
]


# ------------------------------------------------------------
# 20. Convert boundaries to sf format
# ------------------------------------------------------------

kpk_districts_sf <- st_as_sf(
  kpk_districts
)


# ============================================================
# PART VII: CREATE SPATIAL POINT DATA
# ============================================================


# ------------------------------------------------------------
# 21. Convert sampling locations to sf
# ------------------------------------------------------------

mlg_points_sf <- st_as_sf(
  mlg_map_final,
  coords = c(
    "Longitude_round",
    "Latitude_round"
  ),
  crs = 4326
)


# ============================================================
# PART VIII: CREATE GEOGRAPHIC MLG DISTRIBUTION MAP
# ============================================================


# ------------------------------------------------------------
# 22. Create spatial MLG distribution map
# ------------------------------------------------------------

mlg_distribution_map <- ggplot() +
  
  geom_sf(
    data = kpk_districts_sf,
    fill = "grey95",
    color = "grey40",
    linewidth = 0.3
  ) +
  
  geom_sf(
    data = mlg_points_sf,
    aes(
      color = MLG_group,
      size = N_isolates
    ),
    alpha = 0.9
  ) +
  
  scale_size_continuous(
    range = c(3, 10)
  ) +
  
  labs(
    color = "MLG",
    size = "Number of isolates"
  ) +
  
  coord_sf(
    xlim = c(71, 74),
    ylim = c(33.4, 34.7)
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    legend.position = "right",
    
    legend.title = element_text(
      size = 10,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 9
    ),
    
    legend.key.height = unit(
      0.35,
      "cm"
    ),
    
    legend.key.width = unit(
      0.35,
      "cm"
    ),
    
    legend.spacing.y = unit(
      0.15,
      "cm"
    ),
    
    legend.box.spacing = unit(
      0.2,
      "cm"
    ),
    
    axis.title = element_blank()
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(
        size = 4
      )
    ),
    
    size = guide_legend(
      override.aes = list(
        size = c(3, 5, 7, 9)
      )
    )
  )


# ============================================================
# PART IX: CREATE SUPPLEMENTARY TABLES
# ============================================================


# ------------------------------------------------------------
# 23. MLG distribution by district
# ------------------------------------------------------------

MLG_district_table <- mlg_geo %>%
  count(
    MLG,
    District,
    name = "Number_of_isolates"
  )


# ------------------------------------------------------------
# 24. Isolate-level MLG assignments
# ------------------------------------------------------------

MLG_isolate_table <- mlg_geo %>%
  select(
    Isolate,
    District,
    Locality,
    `Longitude (E)`,
    `Latitude (N)`,
    MLG
  )


# ============================================================
# PART X: EXPORT RESULTS
# ============================================================


# ------------------------------------------------------------
# 25. Export MLG distribution by district
# ------------------------------------------------------------

write.csv(
  MLG_district_table,
  file = file.path(
    output_dir,
    "Supplementary_Table_MLG_distribution_by_district.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 26. Export isolate-level MLG assignments
# ------------------------------------------------------------

write.csv(
  MLG_isolate_table,
  file = file.path(
    output_dir,
    "Supplementary_Table_Isolate_MLG_assignment.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 27. Export MLG frequency summary
# ------------------------------------------------------------

write.csv(
  mlg_frequency,
  file = file.path(
    output_dir,
    "Supplementary_Table_MLG_frequency_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 28. Export publication-quality map
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "MLG_geographic_distribution_KPK.tiff"
  ),
  plot = mlg_distribution_map,
  device = "tiff",
  dpi = 600,
  compression = "lzw",
  width = 8,
  height = 6,
  units = "in"
)


# ------------------------------------------------------------
# 29. Report summary
# ------------------------------------------------------------

cat(
  "\nSpatial MLG distribution analysis complete.\n\n"
)

cat(
  "Total isolates:",
  nrow(mlg_geo),
  "\n"
)

cat(
  "Total MLGs:",
  n_distinct(mlg_geo$MLG),
  "\n"
)

cat(
  "Dominant MLGs (>= 5 isolates):",
  length(dominant_MLGs),
  "\n"
)

cat(
  "Dominant genotypes:",
  paste(
    dominant_MLGs,
    collapse = ", "
  ),
  "\n"
)

cat(
  "\nResults exported to:",
  output_dir,
  "\n"
)

print(
  mlg_frequency
)

print(
  mlg_distribution_map
)