# ============================================================
# STRUCTURE Population Structure Visualization
#
# Study:
# Puccinia striiformis f. sp. tritici populations from
# Khyber Pakhtunkhwa, Pakistan
#
# STRUCTURE analysis:
# Bayesian population structure analysis was performed separately
# using STRUCTURE v2.3.4. Replicate runs were aligned and
# summarized using CLUMPAK.
#
# This script does NOT run STRUCTURE or CLUMPAK.
# It imports CLUMPAK-aligned ancestry coefficient files and
# generates STRUCTURE bar plots for K = 2 to K = 7.
#
# Required input files:
# data/structure/
#   K2ClumppIndFile.output
#   K3ClumppIndFile.output
#   K4ClumppIndFile.output
#   K5ClumppIndFile.output
#   K6ClumppIndFile.output
#   K7ClumppIndFile.output
#   Map.xlsx
#
# The final publication figure is provided in the manuscript
# and is not exported by this script.
# ============================================================


# ============================================================
# 1. Load required packages
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)


# ============================================================
# 2. Define input directory
# ============================================================

input_dir <- file.path(
  "data",
  "structure"
)


# ============================================================
# 3. Function to import CLUMPAK-aligned ancestry coefficients
# ============================================================

read_clumpp <- function(file) {
  
  raw <- readLines(file)
  
  ancestry_matrix <- do.call(
    rbind,
    lapply(
      raw,
      function(x) {
        as.numeric(
          strsplit(
            trimws(sub(".*:", "", x)),
            "\\s+"
          )[[1]]
        )
      }
    )
  )
  
  ancestry_df <- as.data.frame(ancestry_matrix)
  
  colnames(ancestry_df) <- paste0(
    "Cluster_",
    seq_len(ncol(ancestry_df))
  )
  
  return(ancestry_df)
}


# ============================================================
# 4. Import CLUMPAK-aligned results for K = 2 to K = 7
# ============================================================

k_values <- 2:7

k_files <- file.path(
  input_dir,
  paste0(
    "K",
    k_values,
    "ClumppIndFile.output"
  )
)

missing_files <- k_files[
  !file.exists(k_files)
]

if (length(missing_files) > 0) {
  
  stop(
    "The following CLUMPAK output files are missing:\n",
    paste(missing_files, collapse = "\n")
  )
}

structure_list <- lapply(
  k_files,
  read_clumpp
)

names(structure_list) <- paste0(
  "K",
  k_values
)


# ============================================================
# 5. Import isolate and population mapping information
# ============================================================

map_file <- file.path(
  input_dir,
  "Map.xlsx"
)

if (!file.exists(map_file)) {
  
  stop(
    "Mapping file not found: ",
    map_file
  )
}

map <- read_excel(map_file)


# ============================================================
# 6. Verify required columns
# ============================================================

required_columns <- c(
  "ID",
  "Population"
)

missing_columns <- required_columns[
  !required_columns %in% names(map)
]

if (length(missing_columns) > 0) {
  
  stop(
    "The following columns are missing from Map.xlsx: ",
    paste(missing_columns, collapse = ", ")
  )
}


# ============================================================
# 7. Define population order
# ============================================================

population_order <- c(
  "Peshawar",
  "Kohat",
  "Charsadda",
  "Mardan",
  "Swabi",
  "Manshera",
  "Buner"
)

map <- map %>%
  mutate(
    Population = factor(
      Population,
      levels = population_order
    )
  ) %>%
  arrange(
    Population,
    ID
  )


# ============================================================
# 8. Validate number of individuals
# ============================================================

n_individuals <- nrow(map)

invalid_k <- names(structure_list)[
  vapply(
    structure_list,
    nrow,
    integer(1)
  ) != n_individuals
]

if (length(invalid_k) > 0) {
  
  stop(
    "The number of individuals in the following CLUMPAK ",
    "output files does not match Map.xlsx: ",
    paste(invalid_k, collapse = ", ")
  )
}


# ============================================================
# 9. Calculate population boundaries
# ============================================================

population_boundaries <- map %>%
  group_by(Population) %>%
  summarise(
    start = min(ID),
    end = max(ID),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(Population)

pop_lines <- population_boundaries$end[
  -nrow(population_boundaries)
]


# ============================================================
# 10. Define ancestry colours
# ============================================================

structure_colors <- c(
  "Cluster_1" = "#1B9E77",
  "Cluster_2" = "#D95F02",
  "Cluster_3" = "#7570B3",
  "Cluster_4" = "#E7298A",
  "Cluster_5" = "#66A61E",
  "Cluster_6" = "#E6AB02",
  "Cluster_7" = "#1F78B4"
)


# ============================================================
# 11. Validate ancestry coefficients
# ============================================================

ancestry_diagnostics <- lapply(
  names(structure_list),
  function(k_name) {
    
    ancestry_df <- structure_list[[k_name]]
    
    data.frame(
      K = k_name,
      Min_value = min(
        as.matrix(ancestry_df),
        na.rm = TRUE
      ),
      Max_value = max(
        as.matrix(ancestry_df),
        na.rm = TRUE
      ),
      Min_row_sum = min(
        rowSums(ancestry_df, na.rm = TRUE)
      ),
      Max_row_sum = max(
        rowSums(ancestry_df, na.rm = TRUE)
      ),
      Missing_values = sum(
        is.na(ancestry_df)
      )
    )
  }
)

ancestry_diagnostics <- do.call(
  rbind,
  ancestry_diagnostics
)

print(ancestry_diagnostics)


# ============================================================
# 12. Create STRUCTURE bar plot function
# ============================================================

make_structure_plot <- function(k) {
  
  ancestry_df <- structure_list[[
    paste0("K", k)
  ]]
  
  
  # Normalize ancestry coefficients so that membership
  # proportions sum to one for each individual
  
  ancestry_df <- ancestry_df %>%
    mutate(
      across(
        everything(),
        ~ .x / rowSums(ancestry_df)
      )
    )
  
  
  # Add individual and population information
  
  plot_df <- ancestry_df %>%
    mutate(
      ID = map$ID,
      Population = map$Population
    ) %>%
    pivot_longer(
      cols = starts_with("Cluster_"),
      names_to = "Cluster",
      values_to = "Ancestry"
    )
  
  
  # Create STRUCTURE bar plot
  
  ggplot(
    plot_df,
    aes(
      x = ID,
      y = Ancestry,
      fill = Cluster
    )
  ) +
    
    geom_bar(
      stat = "identity",
      width = 1
    ) +
    
    geom_vline(
      xintercept = pop_lines + 0.5,
      color = "grey30",
      linewidth = 0.4
    ) +
    
    annotate(
      "text",
      x = -5,
      y = 0.5,
      label = paste0("K = ", k),
      size = 5,
      fontface = "bold",
      hjust = 0.5
    ) +
    
    scale_fill_manual(
      values = structure_colors
    ) +
    
    scale_x_continuous(
      limits = c(
        -10,
        n_individuals
      ),
      expand = c(0, 0)
    ) +
    
    scale_y_continuous(
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    
    theme_classic(
      base_size = 14
    ) +
    
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.6
      ),
      
      plot.margin = margin(
        5,
        5,
        5,
        5
      )
    )
}


# ============================================================
# 13. Generate STRUCTURE plots for K = 2 to K = 7
# ============================================================

structure_plots <- lapply(
  k_values,
  make_structure_plot
)

names(structure_plots) <- paste0(
  "K",
  k_values
)


# ============================================================
# 14. Combine plots vertically
# ============================================================

structure_combined <- wrap_plots(
  structure_plots,
  ncol = 1
)

print(structure_combined)


# ============================================================
# 15. Report completion
# ============================================================

cat(
  "\nSTRUCTURE visualization completed successfully.\n"
)

cat(
  "Individuals:",
  n_individuals,
  "\n"
)

cat(
  "Populations:",
  nlevels(map$Population),
  "\n"
)

cat(
  "STRUCTURE solutions plotted:",
  paste(
    paste0("K = ", k_values),
    collapse = ", "
  ),
  "\n"
)

cat(
  "The final publication figure is provided in the manuscript.\n"
)