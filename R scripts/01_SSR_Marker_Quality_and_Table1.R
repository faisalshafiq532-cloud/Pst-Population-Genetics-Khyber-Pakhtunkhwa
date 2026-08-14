# ============================================================
# SSR Marker Quality Assessment and Generation of Table 1
#
# Study: Puccinia striiformis f. sp. tritici populations from
#        Khyber Pakhtunkhwa, Pakistan
#
# Input:
#   SSR_PK-data.xlsx
#   Sheet: PK
#
# Output:
#   Table_1_SSR_Marker_Quality.csv
#   Table_1_SSR_Marker_Quality.xlsx
#
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(readxl)
library(adegenet)
library(PopGenReport)
library(openxlsx)


# ------------------------------------------------------------
# 2. Import SSR genotype dataset
# ------------------------------------------------------------

input_file <- "SSR_PK-data.xlsx"

ssr <- read_excel(input_file, sheet = "PK")


# ------------------------------------------------------------
# 3. Define SSR loci and validate genotype coding
# ------------------------------------------------------------

# The 17 SSR loci are located in columns 4-20.
locus_cols <- names(ssr)[4:20]

# Check that all non-missing genotypes contain two concatenated
# three-digit allele sizes (e.g., 193195).
genotype_check <- data.frame(
  Locus = locus_cols,
  Missing = sapply(ssr[locus_cols], function(x) sum(is.na(x))),
  Non_6_digit = sapply(ssr[locus_cols], function(x) {
    vals <- as.character(x[!is.na(x)])
    sum(nchar(vals) != 6)
  })
)

# Stop the analysis if invalid genotype entries are detected.
if (sum(genotype_check$Non_6_digit) > 0) {
  stop(
    "Invalid SSR genotype coding detected. ",
    "All non-missing genotypes must contain six digits ",
    "(two three-digit allele sizes)."
  )
}


# ------------------------------------------------------------
# 4. Determine allele-size ranges for each SSR locus
# ------------------------------------------------------------

allele_range_check <- data.frame(
  Locus = locus_cols,
  Min_allele = sapply(ssr[locus_cols], function(x) {
    
    g <- as.character(x)
    g <- g[!is.na(g)]
    
    a1 <- as.numeric(substr(g, 1, 3))
    a2 <- as.numeric(substr(g, 4, 6))
    
    min(c(a1, a2))
  }),
  
  Max_allele = sapply(ssr[locus_cols], function(x) {
    
    g <- as.character(x)
    g <- g[!is.na(g)]
    
    a1 <- as.numeric(substr(g, 1, 3))
    a2 <- as.numeric(substr(g, 4, 6))
    
    max(c(a1, a2))
  })
)


# ------------------------------------------------------------
# 5. Convert concatenated genotypes to allele-pair format
# ------------------------------------------------------------

# Example:
# 193195 -> 193/195

ssr_genotypes <- ssr[locus_cols]

for (locus in locus_cols) {
  
  g <- as.character(ssr_genotypes[[locus]])
  
  allele1 <- substr(g, 1, 3)
  allele2 <- substr(g, 4, 6)
  
  ssr_genotypes[[locus]] <- paste(allele1, allele2, sep = "/")
}


# ------------------------------------------------------------
# 6. Verify genotype integrity after conversion
# ------------------------------------------------------------

genotype_integrity <- data.frame(
  Locus = locus_cols,
  
  Genotypes = sapply(ssr_genotypes, function(x) {
    sum(!is.na(x))
  }),
  
  Invalid_genotypes = sapply(ssr_genotypes, function(x) {
    sum(!is.na(x) & !grepl("^[0-9]{3}/[0-9]{3}$", x))
  })
)

if (sum(genotype_integrity$Invalid_genotypes) > 0) {
  stop(
    "Invalid genotype entries detected after allele-pair conversion."
  )
}


# ------------------------------------------------------------
# 7. Convert SSR data to a genind object
# ------------------------------------------------------------

gen <- df2genind(
  ssr_genotypes,
  sep = "/",
  ploidy = 2,
  ind.names = paste0("PK_", seq_len(nrow(ssr))),
  pop = ssr$Regions,
  type = "codom"
)


# ------------------------------------------------------------
# 8. Calculate observed number of alleles (Na)
# ------------------------------------------------------------

Na <- nAll(gen)


# ------------------------------------------------------------
# 9. Calculate effective number of alleles (Ne)
#
# Ne = 1 / sum(p_i^2)
# ------------------------------------------------------------

# Extract allele-count matrix.
allele_counts <- tab(gen, freq = FALSE)

# Calculate global allele frequencies across all isolates.
allele_freq_global <- colSums(allele_counts) / (2 * nInd(gen))

# Calculate Ne separately for each locus.
Ne <- sapply(levels(locFac(gen)), function(locus) {
  
  p <- allele_freq_global[locFac(gen) == locus]
  p <- p[p > 0]
  
  1 / sum(p^2)
})


# ------------------------------------------------------------
# 10. Calculate observed heterozygosity (Ho)
# ------------------------------------------------------------

Ho <- sapply(locus_cols, function(locus) {
  
  g <- as.character(ssr_genotypes[[locus]])
  valid <- !is.na(g)
  
  alleles <- strsplit(g[valid], "/", fixed = TRUE)
  
  heterozygous <- sapply(
    alleles,
    function(x) x[1] != x[2]
  )
  
  mean(heterozygous)
})


# ------------------------------------------------------------
# 11. Calculate expected heterozygosity (He)
#
# He = 1 - sum(p_i^2)
# ------------------------------------------------------------

He <- sapply(levels(locFac(gen)), function(locus) {
  
  p <- allele_freq_global[locFac(gen) == locus]
  p <- p[p > 0]
  
  1 - sum(p^2)
})


# ------------------------------------------------------------
# 12. Calculate Shannon information index (I)
#
# I = -sum(p_i * ln(p_i))
# ------------------------------------------------------------

Shannon_I <- sapply(levels(locFac(gen)), function(locus) {
  
  p <- allele_freq_global[locFac(gen) == locus]
  p <- p[p > 0]
  
  -sum(p * log(p))
})


# ------------------------------------------------------------
# 13. Calculate polymorphism information content (PIC)
# ------------------------------------------------------------

PIC <- sapply(levels(locFac(gen)), function(locus) {
  
  p <- allele_freq_global[locFac(gen) == locus]
  p <- p[p > 0]
  
  term1 <- sum(p^2)
  
  pairwise_terms <- outer(p^2, p^2)
  term2 <- 2 * sum(pairwise_terms[lower.tri(pairwise_terms)])
  
  1 - term1 - term2
})


# ------------------------------------------------------------
# 14. Estimate null-allele frequencies
#
# Brookfield estimator based on null.all() output.
# ------------------------------------------------------------

null_results <- null.all(gen)

null_brookfield <- as.data.frame(
  t(null_results$null.allele.freq$summary2)
)

null_brookfield$Locus <- toupper(rownames(null_brookfield))

null_brookfield <- null_brookfield[, c(
  "Locus",
  "Observed frequency",
  "Median frequency",
  "2.5th percentile",
  "97.5th percentile"
)]

rownames(null_brookfield) <- NULL


# ------------------------------------------------------------
# 15. Estimate null-allele frequencies within each district
# ------------------------------------------------------------

pop_names <- levels(pop(gen))

null_by_pop <- lapply(pop_names, function(p) {
  
  gen_pop <- gen[pop(gen) == p]
  
  res <- null.all(gen_pop)
  
  estimates <- res$null.allele.freq$summary2[
    "Observed frequency",
  ]
  
  as.numeric(estimates)
})

null_by_pop_table <- as.data.frame(
  do.call(rbind, null_by_pop)
)

rownames(null_by_pop_table) <- pop_names
colnames(null_by_pop_table) <- toupper(locus_cols)


# ------------------------------------------------------------
# 16. Create publication-ready Table 1
# ------------------------------------------------------------

Table1_SSR <- data.frame(
  
  Locus = toupper(locus_cols),
  
  `Allele range (bp)` = paste0(
    allele_range_check$Min_allele,
    "\u2013",
    allele_range_check$Max_allele
  ),
  
  Alleles = as.integer(Na),
  Ne = round(as.numeric(Ne), 2),
  Ho = round(as.numeric(Ho), 3),
  He = round(as.numeric(He), 3),
  Shannon_I = round(as.numeric(Shannon_I), 3),
  PIC = round(as.numeric(PIC), 3),
  `Missing (%)` = 0,
  
  check.names = FALSE
)


# ------------------------------------------------------------
# 17. Add mean row
# ------------------------------------------------------------

Table1_mean <- data.frame(
  
  Locus = "Mean",
  `Allele range (bp)` = "",
  
  Alleles = round(mean(Table1_SSR$Alleles), 2),
  Ne = round(mean(Table1_SSR$Ne), 2),
  Ho = round(mean(Table1_SSR$Ho), 3),
  He = round(mean(Table1_SSR$He), 3),
  Shannon_I = round(mean(Table1_SSR$Shannon_I), 3),
  PIC = round(mean(Table1_SSR$PIC), 3),
  `Missing (%)` = 0,
  
  check.names = FALSE
)

Table1_SSR_final <- rbind(
  Table1_SSR,
  Table1_mean
)


# ------------------------------------------------------------
# 18. Export Table 1
# ------------------------------------------------------------

write.csv(
  Table1_SSR_final,
  "Table_1_SSR_Marker_Quality.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.xlsx(
  Table1_SSR_final,
  "Table_1_SSR_Marker_Quality.xlsx",
  rowNames = FALSE
)