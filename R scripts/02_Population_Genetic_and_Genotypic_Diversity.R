# ============================================================
# Population Genetic and Genotypic Diversity Analysis
#
# Study: Puccinia striiformis f. sp. tritici populations from
#        Khyber Pakhtunkhwa, Pakistan
#
# Input:
#   SSR_PK-data.xlsx
#   Sheet: PK
#
# Main outputs:
#   Population_Genetic_Diversity_Main_Table.csv
#   Population_Genetic_Diversity_Main_Table.xlsx
#
# Supplementary outputs:
#   Population_Genetic_Diversity_Supplementary.xlsx
#
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(readxl)
library(adegenet)
library(poppr)
library(vegan)
library(hierfstat)
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
  
  g <- as.character(ssr_genotypes[[locus]])
  
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

mlg_vector <- poppr::mlg.vector(gen)

mlg_frequency <- sort(
  table(mlg_vector),
  decreasing = TRUE
)

total_MLGs <- length(
  unique(mlg_vector)
)

singleton_MLGs <- sum(
  mlg_frequency == 1
)

repeated_MLGs <- sum(
  mlg_frequency > 1
)


# ------------------------------------------------------------
# 6. Generate the MLG-by-population frequency matrix
# ------------------------------------------------------------

mlg_table_pop <- poppr::mlg.table(
  gen,
  plot = FALSE,
  total = FALSE
)

# Apply the predefined district order.
mlg_table_pop <- mlg_table_pop[
  district_order,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 7. Calculate population-level genotypic diversity
# ------------------------------------------------------------

genotypic_diversity <- poppr::diversity_stats(
  mlg_table_pop,
  H = TRUE,
  G = TRUE,
  lambda = TRUE,
  E5 = TRUE
)


# ------------------------------------------------------------
# 8. Calculate expected MLG richness (eMLG)
#
# Standardized to the smallest population sample size (n = 7).
# ------------------------------------------------------------

eMLG_n7 <- vegan::rarefy(
  mlg_table_pop,
  sample = 7
)

eMLG_n7 <- as.numeric(eMLG_n7)

names(eMLG_n7) <- rownames(
  mlg_table_pop
)


# ------------------------------------------------------------
# 9. Calculate observed MLG richness and MLG/N ratio
# ------------------------------------------------------------

population_N <- rowSums(
  mlg_table_pop
)

population_MLG <- rowSums(
  mlg_table_pop > 0
)

MLG_diversity_summary <- data.frame(
  Population = rownames(mlg_table_pop),
  N = as.numeric(population_N),
  MLGs = as.numeric(population_MLG),
  MLG_N = round(
    population_MLG / population_N,
    3
  ),
  eMLG_n7 = round(
    eMLG_n7,
    3
  ),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 10. Convert genind data to hierfstat format
# ------------------------------------------------------------

hierfstat_gen <- hierfstat::genind2hierfstat(
  gen
)


# ------------------------------------------------------------
# 11. Calculate rarefied allelic richness (Ar)
#
# Standardized to 14 gene copies = 7 diploid isolates.
# ------------------------------------------------------------

allelic_richness_result <- hierfstat::allelic.richness(
  hierfstat_gen,
  min.n = 14,
  diploid = TRUE
)

Ar_matrix <- allelic_richness_result$Ar

# Ensure district order is consistent.
Ar_matrix <- Ar_matrix[
  ,
  district_order,
  drop = FALSE
]

allelic_richness_summary <- data.frame(
  Population = colnames(Ar_matrix),
  Mean_Ar = apply(
    Ar_matrix,
    2,
    mean,
    na.rm = TRUE
  ),
  SD_Ar = apply(
    Ar_matrix,
    2,
    sd,
    na.rm = TRUE
  ),
  stringsAsFactors = FALSE
)

allelic_richness_summary$Mean_Ar <- round(
  allelic_richness_summary$Mean_Ar,
  3
)

allelic_richness_summary$SD_Ar <- round(
  allelic_richness_summary$SD_Ar,
  3
)


# ------------------------------------------------------------
# 12. Identify private alleles
#
# A private allele is observed in one district and absent from
# all other sampled districts.
# ------------------------------------------------------------

private_alleles <- matrix(
  0,
  nrow = nPop(gen),
  ncol = nLoc(gen),
  dimnames = list(
    levels(pop(gen)),
    locNames(gen)
  )
)

for (locus in seq_len(nLoc(gen))) {
  
  locus_alleles <- gen@tab[
    ,
    gen@loc.fac == levels(gen@loc.fac)[locus],
    drop = FALSE
  ]
  
  allele_present <- locus_alleles > 0
  
  for (p in seq_len(nPop(gen))) {
    
    pop_indices <- which(
      pop(gen) == levels(pop(gen))[p]
    )
    
    other_indices <- which(
      pop(gen) != levels(pop(gen))[p]
    )
    
    in_population <- colSums(
      allele_present[
        pop_indices,
        ,
        drop = FALSE
      ]
    ) > 0
    
    outside_population <- colSums(
      allele_present[
        other_indices,
        ,
        drop = FALSE
      ]
    ) > 0
    
    private_alleles[p, locus] <- sum(
      in_population & !outside_population
    )
  }
}

private_allele_totals <- rowSums(
  private_alleles
)

private_allele_totals <- private_allele_totals[
  district_order
]


# ------------------------------------------------------------
# 13. Independently summarize allele occurrence among districts
# ------------------------------------------------------------

private_check <- list()

for (locus in seq_len(nLoc(gen))) {
  
  locus_name <- locNames(gen)[locus]
  
  locus_data <- gen@tab[
    ,
    gen@loc.fac == levels(gen@loc.fac)[locus],
    drop = FALSE
  ]
  
  allele_names <- colnames(locus_data)
  
  district_occurrence <- sapply(
    seq_along(allele_names),
    function(a) {
      
      sum(
        sapply(
          levels(pop(gen)),
          function(p) {
            
            inds <- which(
              pop(gen) == p
            )
            
            sum(
              locus_data[inds, a]
            ) > 0
          }
        )
      )
    }
  )
  
  private_check[[locus_name]] <- data.frame(
    Locus = toupper(locus_name),
    Allele = sub(
      paste0(
        "^",
        locus_name,
        "\\."
      ),
      "",
      allele_names
    ),
    Number_of_Districts = district_occurrence,
    row.names = NULL
  )
}

private_check_table <- do.call(
  rbind,
  private_check
)

rownames(private_check_table) <- NULL

private_candidates <- subset(
  private_check_table,
  Number_of_Districts == 1
)


# ------------------------------------------------------------
# 14. Calculate Nei's gene diversity (Hs)
# ------------------------------------------------------------

basic_stats <- hierfstat::basic.stats(
  hierfstat_gen,
  diploid = TRUE,
  digits = 6
)

Hs_matrix <- basic_stats$Hs

# Ensure district order is consistent.
Hs_matrix <- Hs_matrix[
  ,
  district_order,
  drop = FALSE
]

Hs_summary <- data.frame(
  Population = colnames(Hs_matrix),
  Nei_Hs = apply(
    Hs_matrix,
    2,
    mean,
    na.rm = TRUE
  ),
  SD_Hs = apply(
    Hs_matrix,
    2,
    sd,
    na.rm = TRUE
  ),
  stringsAsFactors = FALSE
)

Hs_summary$Nei_Hs <- round(
  Hs_summary$Nei_Hs,
  4
)

Hs_summary$SD_Hs <- round(
  Hs_summary$SD_Hs,
  4
)


# ------------------------------------------------------------
# 15. Construct the main manuscript diversity table
# ------------------------------------------------------------

diversity_main <- data.frame(
  Population = district_order,
  N = as.numeric(
    population_N[district_order]
  ),
  MLGs = as.numeric(
    population_MLG[district_order]
  ),
  `MLG/N` = round(
    population_MLG[district_order] /
      population_N[district_order],
    3
  ),
  `eMLG (n=7)` = round(
    eMLG_n7[district_order],
    3
  ),
  `Shannon H` = round(
    genotypic_diversity[district_order, "H"],
    4
  ),
  `Stoddart-Taylor G` = round(
    genotypic_diversity[district_order, "G"],
    4
  ),
  E5 = round(
    genotypic_diversity[district_order, "E5"],
    4
  ),
  `Nei's gene diversity (Hs)` = round(
    Hs_summary$Nei_Hs[
      match(
        district_order,
        Hs_summary$Population
      )
    ],
    4
  ),
  check.names = FALSE
)


# ------------------------------------------------------------
# 16. Construct the extended supplementary diversity table
# ------------------------------------------------------------

diversity_supplementary <- data.frame(
  Population = district_order,
  N = as.numeric(
    population_N[district_order]
  ),
  MLGs = as.numeric(
    population_MLG[district_order]
  ),
  `MLG/N` = round(
    population_MLG[district_order] /
      population_N[district_order],
    3
  ),
  `eMLG (n=7)` = round(
    eMLG_n7[district_order],
    3
  ),
  `Simpson lambda` = round(
    genotypic_diversity[district_order, "lambda"],
    4
  ),
  `Mean allelic richness (Ar)` =
    allelic_richness_summary$Mean_Ar[
      match(
        district_order,
        allelic_richness_summary$Population
      )
    ],
  `SD Ar` =
    allelic_richness_summary$SD_Ar[
      match(
        district_order,
        allelic_richness_summary$Population
      )
    ],
  `Nei Hs` =
    Hs_summary$Nei_Hs[
      match(
        district_order,
        Hs_summary$Population
      )
    ],
  `SD Hs` =
    Hs_summary$SD_Hs[
      match(
        district_order,
        Hs_summary$Population
      )
    ],
  `Private alleles` = as.numeric(
    private_allele_totals[district_order]
  ),
  check.names = FALSE
)


# ------------------------------------------------------------
# 17. Prepare locus-level allelic richness matrix
# ------------------------------------------------------------

Ar_locus_table <- data.frame(
  Locus = toupper(
    rownames(Ar_matrix)
  ),
  Ar_matrix,
  check.names = FALSE
)

colnames(Ar_locus_table)[-1] <- district_order


# ------------------------------------------------------------
# 18. Create output directory
# ------------------------------------------------------------

output_dir <- "Population_Genetic_Analysis"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# ------------------------------------------------------------
# 19. Export main manuscript table
# ------------------------------------------------------------

write.csv(
  diversity_main,
  file = file.path(
    output_dir,
    "Population_Genetic_Diversity_Main_Table.csv"
  ),
  row.names = FALSE
)

write.xlsx(
  diversity_main,
  file = file.path(
    output_dir,
    "Population_Genetic_Diversity_Main_Table.xlsx"
  ),
  rowNames = FALSE
)


# ------------------------------------------------------------
# 20. Export supplementary diversity workbook
# ------------------------------------------------------------

wb <- createWorkbook()

addWorksheet(
  wb,
  "Population_Diversity"
)

addWorksheet(
  wb,
  "Locus_Level_Ar"
)

addWorksheet(
  wb,
  "Private_Alleles"
)

# Population-level supplementary statistics
writeData(
  wb,
  sheet = "Population_Diversity",
  x = diversity_supplementary
)

# Locus-level rarefied allelic richness
writeData(
  wb,
  sheet = "Locus_Level_Ar",
  x = Ar_locus_table
)

# Private-allele assessment
writeData(
  wb,
  sheet = "Private_Alleles",
  x = private_check_table
)

header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center"
)

for (sheet_name in names(wb)) {
  
  addStyle(
    wb,
    sheet = sheet_name,
    style = header_style,
    rows = 1,
    cols = 1:ncol(
      readWorkbook(
        wb,
        sheet = sheet_name
      )
    ),
    gridExpand = TRUE
  )
  
  setColWidths(
    wb,
    sheet = sheet_name,
    cols = 1:50,
    widths = "auto"
  )
  
  freezePane(
    wb,
    sheet = sheet_name,
    firstRow = TRUE
  )
}

saveWorkbook(
  wb,
  file = file.path(
    output_dir,
    "Population_Genetic_Diversity_Supplementary.xlsx"
  ),
  overwrite = TRUE
)


# ------------------------------------------------------------
# 21. Report key results
# ------------------------------------------------------------

cat(
  "\nPopulation genetic and genotypic diversity analysis complete.\n\n"
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
  "Singleton MLGs:",
  singleton_MLGs,
  "\n"
)

cat(
  "Repeated MLGs:",
  repeated_MLGs,
  "\n"
)

cat(
  "Private alleles detected:",
  nrow(private_candidates),
  "\n"
)

print(
  diversity_main
)