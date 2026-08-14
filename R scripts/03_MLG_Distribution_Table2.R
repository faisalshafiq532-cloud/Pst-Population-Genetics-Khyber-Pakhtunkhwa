# ============================================================
# Multilocus Genotype Distribution and Generation of Table 2
#
# Study: Puccinia striiformis f. sp. tritici populations from
#        Khyber Pakhtunkhwa, Pakistan
#
# Input:
#   SSR_PK-data.xlsx
#   Sheet: PK
#
# Output:
#   Table_2_MLG_Distribution_by_Population.xlsx
#
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(readxl)
library(dplyr)
library(writexl)


# ------------------------------------------------------------
# 2. Import SSR genotype dataset
# ------------------------------------------------------------

input_file <- "SSR_PK-data.xlsx"

ssr_data <- read_excel(
  input_file,
  sheet = "PK"
)

# The 17 SSR loci are located in columns 4-20.
marker_cols <- names(ssr_data)[4:20]


# ------------------------------------------------------------
# 3. Define the manuscript population order
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
# 4. Create multilocus genotype (MLG) assignments
#
# Each unique allelic profile across the 17 SSR loci is assigned
# a unique MLG identifier.
# ------------------------------------------------------------

ssr_data <- ssr_data %>%
  mutate(
    MLG_profile = apply(
      select(., all_of(marker_cols)),
      1,
      paste,
      collapse = "_"
    ),
    MLG_number = match(
      MLG_profile,
      unique(MLG_profile)
    ),
    MLG = paste0(
      "MLG.",
      MLG_number
    )
  )


# ------------------------------------------------------------
# 5. Create isolate-level MLG dataset
# ------------------------------------------------------------

mlg_clean <- ssr_data %>%
  mutate(
    Isolate = paste0(
      "PK_",
      row_number()
    ),
    District = Regions
  ) %>%
  select(
    Isolate,
    District,
    MLG
  )


# ------------------------------------------------------------
# 6. Calculate MLG frequencies within each district
# ------------------------------------------------------------

mlg_district_counts <- mlg_clean %>%
  count(
    District,
    MLG,
    name = "Frequency"
  ) %>%
  mutate(
    MLG_number = as.numeric(
      sub(
        "MLG\\.",
        "",
        MLG
      )
    )
  ) %>%
  arrange(
    District,
    desc(Frequency),
    MLG_number
  )


# ------------------------------------------------------------
# 7. Calculate district-level sample size and MLG richness
# ------------------------------------------------------------

district_mlg_summary <- mlg_clean %>%
  group_by(District) %>%
  summarise(
    No_of_isolates = n(),
    No_of_MLGs = n_distinct(MLG),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 8. Classify MLGs by their geographic distribution
#
# S = Shared among two or more districts
# R = Restricted to one district
# ------------------------------------------------------------

mlg_sharing_status <- mlg_clean %>%
  group_by(MLG) %>%
  summarise(
    No_of_districts = n_distinct(District),
    Status = if_else(
      No_of_districts >= 2,
      "S",
      "R"
    ),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 9. Create complete MLG list for each district
#
# Format:
# MLG.1 (frequency; status)
# ------------------------------------------------------------

all_mlg_by_district <- mlg_district_counts %>%
  left_join(
    mlg_sharing_status %>%
      select(
        MLG,
        Status
      ),
    by = "MLG"
  ) %>%
  mutate(
    MLG_entry = paste0(
      MLG,
      " (",
      Frequency,
      "; ",
      Status,
      ")"
    )
  ) %>%
  arrange(
    District,
    desc(Frequency),
    MLG_number
  ) %>%
  group_by(District) %>%
  summarise(
    `MLGs present (frequency; distribution status)` =
      paste(
        MLG_entry,
        collapse = ", "
      ),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 10. Construct manuscript Table 2
# ------------------------------------------------------------

MLG_distribution_final <- district_mlg_summary %>%
  left_join(
    all_mlg_by_district,
    by = "District"
  ) %>%
  mutate(
    District = factor(
      District,
      levels = district_order
    )
  ) %>%
  arrange(District) %>%
  rename(
    Population = District,
    `No. of isolates` = No_of_isolates,
    `No. of MLGs` = No_of_MLGs
  )


# ------------------------------------------------------------
# 11. Add overall total row
#
# Total MLG richness represents the number of unique MLGs
# across the complete dataset, not the sum across populations.
# ------------------------------------------------------------

total_row <- tibble(
  Population = "Total",
  `No. of isolates` = nrow(mlg_clean),
  `No. of MLGs` = n_distinct(mlg_clean$MLG),
  `MLGs present (frequency; distribution status)` = "—"
)

MLG_distribution_final <- bind_rows(
  MLG_distribution_final,
  total_row
)


# ------------------------------------------------------------
# 12. Integrity checks
# ------------------------------------------------------------

# Confirm the total sample size.
stopifnot(
  sum(
    MLG_distribution_final$`No. of isolates`[
      MLG_distribution_final$Population != "Total"
    ]
  ) == nrow(mlg_clean)
)

# Confirm that MLG frequencies sum to the number of isolates
# within every district.
frequency_check <- mlg_district_counts %>%
  group_by(District) %>%
  summarise(
    Sum_of_MLG_frequencies = sum(Frequency),
    .groups = "drop"
  ) %>%
  left_join(
    district_mlg_summary,
    by = "District"
  )

stopifnot(
  all(
    frequency_check$Sum_of_MLG_frequencies ==
      frequency_check$No_of_isolates
  )
)


# ------------------------------------------------------------
# 13. Create output directory
# ------------------------------------------------------------

output_dir <- "MLG_Analysis"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# ------------------------------------------------------------
# 14. Export manuscript Table 2
# ------------------------------------------------------------

write_xlsx(
  MLG_distribution_final,
  path = file.path(
    output_dir,
    "Table_2_MLG_Distribution_by_Population.xlsx"
  )
)


# ------------------------------------------------------------
# 15. Report key results
# ------------------------------------------------------------

cat(
  "\nMLG distribution analysis complete.\n\n"
)

cat(
  "Total isolates:",
  nrow(mlg_clean),
  "\n"
)

cat(
  "Total unique MLGs:",
  n_distinct(mlg_clean$MLG),
  "\n"
)

cat(
  "Shared MLGs:",
  sum(mlg_sharing_status$Status == "S"),
  "\n"
)

cat(
  "District-restricted MLGs:",
  sum(mlg_sharing_status$Status == "R"),
  "\n"
)

print(
  MLG_distribution_final,
  width = Inf
)