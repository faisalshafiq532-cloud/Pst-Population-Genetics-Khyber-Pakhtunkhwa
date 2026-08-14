# ============================================================
# AMOVA and Pairwise Weir-Cockerham FST Analysis
#
# Study: Puccinia striiformis f. sp. tritici populations from
#        Khyber Pakhtunkhwa, Pakistan
#
# Input:
#   SSR_PK-data.xlsx
#   Sheet: PK
#
# Analyses:
#   1. AMOVA based on Bruvo genetic distances
#   2. Pairwise Weir-Cockerham FST
#   3. Bootstrap 95% confidence intervals for pairwise FST
#
# Outputs:
#   AMOVA_and_Pairwise_FST/
#     AMOVA_Table.csv
#     Pairwise_FST_with_95CI.csv
#     Supplementary_AMOVA_and_Pairwise_FST.xlsx
#
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(readxl)
library(adegenet)
library(poppr)
library(pegas)
library(hierfstat)
library(openxlsx)


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
# Example:
# 193195 -> 193/195
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
# 5. Define manuscript population order
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
# 6. Define SSR repeat lengths
#
# Values must correspond exactly to the order of loci in gen.
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
# PART I: AMOVA
# ============================================================


# ------------------------------------------------------------
# 7. Calculate Bruvo genetic distance
# ------------------------------------------------------------

bruvo_dist <- poppr::bruvo.dist(
  gen,
  replen = repeat_lengths,
  add = TRUE
)


# ------------------------------------------------------------
# 8. Prepare district grouping
# ------------------------------------------------------------

amova_data <- data.frame(
  District = factor(
    pop(gen),
    levels = district_order
  )
)

stopifnot(
  nrow(amova_data) ==
    attr(bruvo_dist, "Size")
)


# ------------------------------------------------------------
# 9. Run AMOVA
#
# One-level hierarchical model:
# genetic variation among districts and within districts.
#
# Permutation test: 999 permutations.
# ------------------------------------------------------------

set.seed(20260809)

amova_pst <- pegas::amova(
  bruvo_dist ~ District,
  data = amova_data,
  nperm = 999,
  is.squared = FALSE
)


# ------------------------------------------------------------
# 10. Extract AMOVA results
# ------------------------------------------------------------

# Display the result structure once if needed:
# str(amova_pst)

amova_results <- as.data.frame(
  amova_pst$tab
)

amova_results$Source <- rownames(
  amova_results
)

rownames(amova_results) <- NULL

# Move Source to the first column.
amova_results <- amova_results[
  ,
  c(
    "Source",
    setdiff(
      names(amova_results),
      "Source"
    )
  ),
  drop = FALSE
]


# ============================================================
# PART II: PAIRWISE WEIR-COCKERHAM FST
# ============================================================


# ------------------------------------------------------------
# 11. Convert genind object to hierfstat format
# ------------------------------------------------------------

hierfstat_gen <- hierfstat::genind2hierfstat(
  gen
)


# ------------------------------------------------------------
# 12. Calculate pairwise Weir-Cockerham FST
# ------------------------------------------------------------

pairwise_fst <- hierfstat::pairwise.WCfst(
  hierfstat_gen,
  diploid = TRUE
)

# Apply manuscript population order.
population_order <- district_order[
  district_order %in% rownames(pairwise_fst)
]

pairwise_fst <- pairwise_fst[
  population_order,
  population_order,
  drop = FALSE
]


# ------------------------------------------------------------
# 13. Bootstrap 95% confidence intervals for pairwise FST
#
# Bootstrap resampling across SSR loci.
# ------------------------------------------------------------

set.seed(20260809)

boot_fst <- hierfstat::boot.ppfst(
  dat = hierfstat_gen,
  nboot = 2000,
  quant = c(0.025, 0.975),
  diploid = TRUE
)

boot_lower <- boot_fst$ll[
  population_order,
  population_order,
  drop = FALSE
]

boot_upper <- boot_fst$ul[
  population_order,
  population_order,
  drop = FALSE
]


# ------------------------------------------------------------
# 14. Create pairwise FST publication table
# ------------------------------------------------------------

population_pairs <- combn(
  population_order,
  2,
  simplify = FALSE
)

pairwise_fst_table <- do.call(
  rbind,
  lapply(
    population_pairs,
    function(pair) {
      
      data.frame(
        Population_1 = pair[1],
        Population_2 = pair[2],
        WC_FST = pairwise_fst[
          pair[1],
          pair[2]
        ],
        Lower_95_CI = boot_lower[
          pair[1],
          pair[2]
        ],
        Upper_95_CI = boot_upper[
          pair[1],
          pair[2]
        ],
        stringsAsFactors = FALSE
      )
    }
  )
)

# Round values for publication output.
pairwise_fst_table$WC_FST <- round(
  pairwise_fst_table$WC_FST,
  4
)

pairwise_fst_table$Lower_95_CI <- round(
  pairwise_fst_table$Lower_95_CI,
  4
)

pairwise_fst_table$Upper_95_CI <- round(
  pairwise_fst_table$Upper_95_CI,
  4
)


# ============================================================
# PART III: EXPORT RESULTS
# ============================================================


# ------------------------------------------------------------
# 15. Create output directory
# ------------------------------------------------------------

output_dir <- "AMOVA_and_Pairwise_FST"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# ------------------------------------------------------------
# 16. Export CSV files
# ------------------------------------------------------------

write.csv(
  amova_results,
  file = file.path(
    output_dir,
    "AMOVA_Table.csv"
  ),
  row.names = FALSE
)

write.csv(
  pairwise_fst_table,
  file = file.path(
    output_dir,
    "Pairwise_FST_with_95CI.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. Create supplementary Excel workbook
# ------------------------------------------------------------

wb <- createWorkbook()


# AMOVA worksheet
addWorksheet(
  wb,
  "AMOVA"
)

writeData(
  wb,
  sheet = "AMOVA",
  x = amova_results
)


# Pairwise FST worksheet
addWorksheet(
  wb,
  "Pairwise_FST"
)

writeData(
  wb,
  sheet = "Pairwise_FST",
  x = pairwise_fst_table
)


# ------------------------------------------------------------
# 18. Format workbook
# ------------------------------------------------------------

header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  wrapText = TRUE
)

for (sheet_name in c(
  "AMOVA",
  "Pairwise_FST"
)) {
  
  n_columns <- ncol(
    readWorkbook(
      wb,
      sheet = sheet_name
    )
  )
  
  addStyle(
    wb,
    sheet = sheet_name,
    style = header_style,
    rows = 1,
    cols = seq_len(n_columns),
    gridExpand = TRUE
  )
  
  freezePane(
    wb,
    sheet = sheet_name,
    firstRow = TRUE
  )
  
  setColWidths(
    wb,
    sheet = sheet_name,
    cols = seq_len(n_columns),
    widths = "auto"
  )
}


# ------------------------------------------------------------
# 19. Save supplementary workbook
# ------------------------------------------------------------

saveWorkbook(
  wb,
  file = file.path(
    output_dir,
    "Supplementary_AMOVA_and_Pairwise_FST.xlsx"
  ),
  overwrite = TRUE
)


# ------------------------------------------------------------
# 20. Report key results
# ------------------------------------------------------------

cat(
  "\nAMOVA and pairwise FST analyses complete.\n\n"
)

cat(
  "Number of isolates:",
  nInd(gen),
  "\n"
)

cat(
  "Number of populations:",
  length(population_order),
  "\n"
)

cat(
  "Number of pairwise FST comparisons:",
  nrow(pairwise_fst_table),
  "\n"
)

cat(
  "AMOVA permutations:",
  999,
  "\n"
)

cat(
  "FST bootstrap replicates:",
  2000,
  "\n\n"
)

print(amova_results)

print(pairwise_fst_table)