# ============================================================
# Multilocus Genotype Distribution and Sharing Analysis
#
# Study: Puccinia striiformis f. sp. tritici populations from
#        Khyber Pakhtunkhwa, Pakistan
#
# Input:
#   SSR_PK-data.xlsx
#   Sheet: PK
#
# Output:
#   Supplementary_MLG_Analysis_Tables.xlsx
#
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(readxl)
library(adegenet)
library(poppr)
library(openxlsx)


# ------------------------------------------------------------
# 2. Import and prepare the SSR dataset
# ------------------------------------------------------------

input_file <- "SSR_PK-data.xlsx"

ssr <- read_excel(
  input_file,
  sheet = "PK"
)

# The 17 SSR loci are located in columns 4-20.
locus_cols <- names(ssr)[4:20]

# Convert concatenated genotypes to allele-pair format.
# Example: 193195 -> 193/195
ssr_genotypes <- ssr[locus_cols]

for (locus in locus_cols) {
  
  g <- as.character(
    ssr_genotypes[[locus]]
  )
  
  allele1 <- substr(g, 1, 3)
  allele2 <- substr(g, 4, 6)
  
  ssr_genotypes[[locus]] <- paste(
    allele1,
    allele2,
    sep = "/"
  )
}


# ------------------------------------------------------------
# 3. Create a genind object
# ------------------------------------------------------------

gen <- df2genind(
  ssr_genotypes,
  sep = "/",
  ploidy = 2,
  ind.names = paste0(
    "PK_",
    seq_len(nrow(ssr))
  ),
  pop = ssr$Regions,
  type = "codom"
)


# ------------------------------------------------------------
# 4. Define district order
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
# 5. Identify multilocus genotypes
# ------------------------------------------------------------

mlg_vector <- poppr::mlg.vector(
  gen
)

mlg_frequency <- sort(
  table(mlg_vector),
  decreasing = TRUE
)


# ------------------------------------------------------------
# 6. Generate the MLG-by-population frequency matrix
# ------------------------------------------------------------

mlg_table_pop <- poppr::mlg.table(
  gen,
  plot = FALSE,
  total = FALSE
)

mlg_table_pop <- mlg_table_pop[
  district_order,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 7. Classify MLGs by geographic distribution
# ------------------------------------------------------------

# Number of districts in which each MLG occurs.
mlg_n_districts <- colSums(
  mlg_table_pop > 0
)

# Total number of isolates represented by each MLG.
mlg_n_isolates <- colSums(
  mlg_table_pop
)

# Classify MLGs.
mlg_class <- ifelse(
  mlg_n_districts >= 2,
  "Shared across districts",
  "District-restricted"
)

MLG_crosspop_summary <- data.frame(
  MLG = colnames(mlg_table_pop),
  Frequency = as.numeric(
    mlg_n_isolates
  ),
  Number_of_Districts = as.numeric(
    mlg_n_districts
  ),
  Classification = mlg_class,
  stringsAsFactors = FALSE
)

# Sort from the most geographically widespread MLGs to the
# least widespread, then by total frequency.
MLG_crosspop_summary <- MLG_crosspop_summary[
  order(
    -MLG_crosspop_summary$Number_of_Districts,
    -MLG_crosspop_summary$Frequency
  ),
  ,
  drop = FALSE
]

rownames(MLG_crosspop_summary) <- NULL


# ------------------------------------------------------------
# 8. Calculate headline MLG sharing statistics
# ------------------------------------------------------------

total_MLGs <- nrow(
  MLG_crosspop_summary
)

shared_MLGs <- sum(
  MLG_crosspop_summary$Number_of_Districts >= 2
)

restricted_MLGs <- sum(
  MLG_crosspop_summary$Number_of_Districts == 1
)

shared_isolates <- sum(
  MLG_crosspop_summary$Frequency[
    MLG_crosspop_summary$Number_of_Districts >= 2
  ]
)

restricted_isolates <- sum(
  MLG_crosspop_summary$Frequency[
    MLG_crosspop_summary$Number_of_Districts == 1
  ]
)

MLG_sharing_summary <- data.frame(
  Metric = c(
    "Total MLGs",
    "Shared MLGs",
    "District-restricted MLGs",
    "Percentage of MLGs shared",
    "Percentage of MLGs district-restricted",
    "Isolates belonging to shared MLGs",
    "Isolates belonging to district-restricted MLGs",
    "Percentage of isolates in shared MLGs",
    "Percentage of isolates in district-restricted MLGs"
  ),
  
  Value = c(
    total_MLGs,
    shared_MLGs,
    restricted_MLGs,
    round(
      100 * shared_MLGs / total_MLGs,
      1
    ),
    round(
      100 * restricted_MLGs / total_MLGs,
      1
    ),
    shared_isolates,
    restricted_isolates,
    round(
      100 * shared_isolates / nInd(gen),
      1
    ),
    round(
      100 * restricted_isolates / nInd(gen),
      1
    )
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 9. Create complete MLG × district distribution table
# ------------------------------------------------------------

MLG_distribution <- as.data.frame(
  mlg_table_pop
)

MLG_distribution$Population <- rownames(
  MLG_distribution
)

MLG_distribution <- MLG_distribution[
  ,
  c(
    "Population",
    setdiff(
      names(MLG_distribution),
      "Population"
    )
  ),
  drop = FALSE
]

rownames(MLG_distribution) <- NULL


# ------------------------------------------------------------
# 10. Create table of MLGs shared among districts
# ------------------------------------------------------------

shared_mlg_names <- MLG_crosspop_summary$MLG[
  MLG_crosspop_summary$Number_of_Districts >= 2
]

MLG_shared_distribution <- data.frame(
  MLG = shared_mlg_names,
  
  Frequency =
    MLG_crosspop_summary$Frequency[
      MLG_crosspop_summary$Number_of_Districts >= 2
    ],
  
  Number_of_Districts =
    MLG_crosspop_summary$Number_of_Districts[
      MLG_crosspop_summary$Number_of_Districts >= 2
    ],
  
  Distribution = NA_character_,
  
  stringsAsFactors = FALSE
)

for (i in seq_len(
  nrow(MLG_shared_distribution)
)) {
  
  this_mlg <- MLG_shared_distribution$MLG[i]
  
  present_districts <- rownames(
    mlg_table_pop
  )[
    mlg_table_pop[, this_mlg] > 0
  ]
  
  MLG_shared_distribution$Distribution[i] <- paste(
    present_districts,
    collapse = "; "
  )
}

MLG_shared_distribution <- MLG_shared_distribution[
  order(
    -MLG_shared_distribution$Number_of_Districts,
    -MLG_shared_distribution$Frequency
  ),
  ,
  drop = FALSE
]

rownames(MLG_shared_distribution) <- NULL


# ------------------------------------------------------------
# 11. Calculate district-level contribution of shared MLGs
# ------------------------------------------------------------

shared_by_pop <- rowSums(
  mlg_table_pop[
    ,
    shared_mlg_names,
    drop = FALSE
  ]
)

population_sizes <- rowSums(
  mlg_table_pop
)

MLG_sharing_by_population <- data.frame(
  Population = rownames(
    mlg_table_pop
  ),
  
  N = as.numeric(
    population_sizes
  ),
  
  Shared_MLG_Isolates = as.numeric(
    shared_by_pop
  ),
  
  Percentage_Shared = round(
    100 * shared_by_pop /
      population_sizes,
    1
  ),
  
  stringsAsFactors = FALSE
)

rownames(MLG_sharing_by_population) <- NULL


# ------------------------------------------------------------
# 12. Create output directory
# ------------------------------------------------------------

output_dir <- "MLG_Analysis"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# ------------------------------------------------------------
# 13. Export headline MLG summary
# ------------------------------------------------------------

write.csv(
  MLG_crosspop_summary,
  file = file.path(
    output_dir,
    "MLG_Cross_Population_Summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  MLG_sharing_summary,
  file = file.path(
    output_dir,
    "MLG_Sharing_Headline_Statistics.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. Create supplementary MLG workbook
# ------------------------------------------------------------

wb <- createWorkbook()

addWorksheet(
  wb,
  "MLG_Distribution"
)

addWorksheet(
  wb,
  "Shared_MLGs"
)

addWorksheet(
  wb,
  "District_MLG_Sharing"
)

addWorksheet(
  wb,
  "MLG_Summary"
)


# ------------------------------------------------------------
# 15. Define workbook styles
# ------------------------------------------------------------

title_style <- createStyle(
  textDecoration = "bold",
  wrapText = TRUE,
  valign = "top"
)

header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom"
)


# ------------------------------------------------------------
# 16. Write complete MLG distribution table
# ------------------------------------------------------------

writeData(
  wb,
  sheet = "MLG_Distribution",
  x = MLG_distribution,
  startRow = 2
)

writeData(
  wb,
  sheet = "MLG_Distribution",
  x = paste(
    "Distribution and frequency of multilocus genotypes",
    "(MLGs) among Pst isolates across seven districts of",
    "Khyber Pakhtunkhwa, Pakistan."
  ),
  startRow = 1,
  startCol = 1
)

addStyle(
  wb,
  sheet = "MLG_Distribution",
  style = title_style,
  rows = 1,
  cols = 1
)

addStyle(
  wb,
  sheet = "MLG_Distribution",
  style = header_style,
  rows = 2,
  cols = 1:ncol(MLG_distribution),
  gridExpand = TRUE
)

freezePane(
  wb,
  sheet = "MLG_Distribution",
  firstRow = TRUE,
  firstCol = TRUE
)


# ------------------------------------------------------------
# 17. Write shared MLG distribution table
# ------------------------------------------------------------

writeData(
  wb,
  sheet = "Shared_MLGs",
  x = MLG_shared_distribution
)

addStyle(
  wb,
  sheet = "Shared_MLGs",
  style = header_style,
  rows = 1,
  cols = 1:ncol(MLG_shared_distribution),
  gridExpand = TRUE
)

freezePane(
  wb,
  sheet = "Shared_MLGs",
  firstRow = TRUE
)


# ------------------------------------------------------------
# 18. Write district-level MLG sharing table
# ------------------------------------------------------------

writeData(
  wb,
  sheet = "District_MLG_Sharing",
  x = MLG_sharing_by_population
)

addStyle(
  wb,
  sheet = "District_MLG_Sharing",
  style = header_style,
  rows = 1,
  cols = 1:ncol(
    MLG_sharing_by_population
  ),
  gridExpand = TRUE
)

freezePane(
  wb,
  sheet = "District_MLG_Sharing",
  firstRow = TRUE
)


# ------------------------------------------------------------
# 19. Write overall MLG summary
# ------------------------------------------------------------

writeData(
  wb,
  sheet = "MLG_Summary",
  x = MLG_crosspop_summary,
  startRow = 1
)

writeData(
  wb,
  sheet = "MLG_Summary",
  x = MLG_sharing_summary,
  startRow = nrow(
    MLG_crosspop_summary
  ) + 4
)

addStyle(
  wb,
  sheet = "MLG_Summary",
  style = header_style,
  rows = 1,
  cols = 1:ncol(
    MLG_crosspop_summary
  ),
  gridExpand = TRUE
)


# ------------------------------------------------------------
# 20. Adjust workbook formatting
# ------------------------------------------------------------

for (sheet_name in names(wb)) {
  
  setColWidths(
    wb,
    sheet = sheet_name,
    cols = 1:50,
    widths = "auto"
  )
}


# ------------------------------------------------------------
# 21. Save supplementary workbook
# ------------------------------------------------------------

saveWorkbook(
  wb,
  file = file.path(
    output_dir,
    "Supplementary_MLG_Analysis_Tables.xlsx"
  ),
  overwrite = TRUE
)


# ------------------------------------------------------------
# 22. Report key results
# ------------------------------------------------------------

cat(
  "\nMLG distribution and sharing analysis complete.\n\n"
)

cat(
  "Total isolates:",
  nInd(gen),
  "\n"
)

cat(
  "Total MLGs:",
  total_MLGs,
  "\n"
)

cat(
  "Shared MLGs:",
  shared_MLGs,
  "\n"
)

cat(
  "District-restricted MLGs:",
  restricted_MLGs,
  "\n"
)

cat(
  "Percentage of MLGs shared:",
  round(
    100 * shared_MLGs / total_MLGs,
    1
  ),
  "%\n"
)

print(
  MLG_sharing_by_population
)