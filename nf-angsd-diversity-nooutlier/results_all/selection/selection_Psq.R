# =========================================================
# Load packages
# =========================================================
library(data.table)
library(dplyr)
library(ACER)

# =========================================================
# Set working directory
# =========================================================
setwd("/archive/carpenterlab/pire/pire_pseudanthias_squamipinnis_lcwgs/nf-angsd-diversity-nooutlier/results_all/angsd_pop/")

# =========================================================
# READ FILES
# =========================================================

agal_mafs <- fread(
  "PsqAGal_historic.mafs.gz",
  header = TRUE
)

cver_mafs <- fread(
  "PsqCVer_modern.mafs.gz",
  header = TRUE
)

# =========================================================
# FIX MAJOR / REF
# =========================================================

fix_major_knownEM <- function(df) {
  
  df %>%
    mutate(
      
      knownEM = if_else(
        major != ref,
        1 - knownEM,
        knownEM
      ),
      
      major = if_else(
        major != ref,
        ref,
        major
      )
    )
}

agal_mafs <- fix_major_knownEM(agal_mafs)
cver_mafs <- fix_major_knownEM(cver_mafs)

# =========================================================
# MAF FILTER: Minimum allele count = 5 / Total individuals = 62
# 5 / 62 = 0.081
# =========================================================

maf_filter <- 0.081

# =========================================================
# KEEP SHARED SNPs + APPLY MAF FILTER
# =========================================================

filtered <- inner_join(
  agal_mafs,
  cver_mafs,
  by = c("chromo", "position"),
  suffix = c("_agal", "_cver")
) %>%
  filter(
    knownEM_agal > maf_filter |
      knownEM_cver > maf_filter
  )

cat(
  "SNPs after MAF filtering:",
  nrow(filtered),
  "\n"
)
#shared SNPs after MAF filtering: 956,656
# =========================================================
# CREATE FREQUENCY MATRIX
# =========================================================

freq_mat <- data.frame(
  agal = filtered$knownEM_agal,
  cver = filtered$knownEM_cver
)

# =========================================================
# CREATE COVERAGE MATRIX
# =========================================================

cov_mat <- data.frame(
  agal = filtered$nInd_agal,
  cver = filtered$nInd_cver
)

# =========================================================
# PARAMETERS
# =========================================================

Ne <- c(1452, 1452)

# generations between historical & modern
gen <- c(0, 63)

# =========================================================
# RUN ADAPTED CHI-SQUARE TEST
# =========================================================

pval <- adapted.chisq.test(
  freq     = freq_mat,
  coverage = cov_mat,
  Ne       = Ne,
  gen      = gen
)

# =========================================================
# CONVERT TO DATAFRAME
# =========================================================

pval <- as.data.frame(pval)

# =========================================================
# FDR CORRECTION
# =========================================================

pval$fdr <- p.adjust(
  pval$pval,
  method = "fdr"
)

# =========================================================
# COMBINE RESULTS
# =========================================================

freq <- cbind(
  pval,
  freq_mat
)

freq$chromo <- filtered$chromo
freq$position <- filtered$position

# =========================================================
# SNPs UNDER SELECTION
# =========================================================

freq_fdr <- subset(
  freq,
  fdr < 0.05
)

# =========================================================
# ALLELE FREQUENCY CHANGE
# =========================================================

freq_fdr$delta <-
  freq_fdr$cver -
  freq_fdr$agal

# =========================================================
# SUMMARY STATISTICS
# =========================================================

cat(
  "Total SNPs tested:",
  nrow(freq),
  "\n"
)
#Total SNPs tested: 956,656

cat(
  "SNPs under selection:",
  nrow(freq_fdr),
  "\n"
)
#SNPs under selection: 12,196

cat(
  "Percent selected:",
  round(
    nrow(freq_fdr) /
      nrow(freq) * 100,
    4
  ),
  "%\n"
)
#Percent selected: 1.2749 %
# =========================================================
# DIRECTION OF CHANGE
# =========================================================

table(
  sign(freq_fdr$delta)
)
#  -1    1 
# 1575 10621 
# =========================================================
# POSITIVE VS NEGATIVE CHANGE
# =========================================================

cat(
  "Positive allele frequency shifts:",
  sum(freq_fdr$delta > 0),
  "\n"
)
#Positive allele frequency shifts: 10,621

cat(
  "Negative allele frequency shifts:",
  sum(freq_fdr$delta < 0),
  "\n"
)
#Negative allele frequency shifts: 1,575 
# =========================================================
# LARGE FREQUENCY SHIFTS
# =========================================================

large_shift <- freq_fdr %>%
  filter(
    abs(delta) > 0.5
  )

cat(
  "SNPs with |delta| > 0.5:",
  nrow(large_shift),
  "\n"
)
#SNPs with |delta| > 0.5: 492  
# -----------------------------
# Neutral SNP extraction
# -----------------------------
# Create chrpos IDs
freq$chrpos <- paste(freq$chromo, freq$position, sep = "_")
freq_fdr$chrpos <- paste(freq_fdr$chromo, freq_fdr$position, sep = "_")

neutral_snps <- freq[!freq$chrpos %in% freq_fdr$chrpos, ]

cat("Neutral SNPs:", nrow(neutral_snps), "\n")
#Globally Neutral SNPs: 944,460

# Export neutral SNPs
write.table(
  neutral_snps[, c("chromo","position")],
  file = "snp_neutral.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

cat("Neutral SNP list saved.\n")


