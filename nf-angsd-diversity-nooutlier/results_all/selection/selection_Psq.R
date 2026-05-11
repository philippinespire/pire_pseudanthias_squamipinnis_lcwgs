# -----------------------------
# Load packages
# -----------------------------
library(data.table)
library(dplyr)
library(ACER)

# -----------------------------
# Set working directory
# -----------------------------
setwd("/archive/carpenterlab/pire/pire_pseudanthias_squamipinnis_lcwgs/nf-angsd-diversity-nooutlier/results/angsd_pop/")

# -----------------------------
# Load data
# -----------------------------
agal_mafs <- fread("PsqAGal_historic.mafs.gz", header=TRUE)
cver_mafs <- fread("PsqCVer_modern.mafs.gz", header=TRUE)

# -----------------------------
# Fix major/ref and knownEM
# -----------------------------
fix_major_knownEM <- function(df) {
  df %>%
    mutate(
      knownEM = if_else(major != ref, 1 - knownEM, knownEM),
      major = if_else(major != ref, ref, major)
    )
}

agal_mafs <- fix_major_knownEM(agal_mafs)
cver_mafs <- fix_major_knownEM(cver_mafs)

# -----------------------------
# STEP 1: Keep only shared SNPs
# -----------------------------
shared_snps <- inner_join(
  agal_mafs[, .(chromo, position)],
  cver_mafs[, .(chromo, position)],
  by = c("chromo", "position")
)

# Filter datasets
agal_mafs_filt <- agal_mafs %>%
  semi_join(shared_snps, by = c("chromo", "position"))

cver_mafs_filt <- cver_mafs %>%
  semi_join(shared_snps, by = c("chromo", "position"))

# -----------------------------
# STEP 2: Ensure identical order
# -----------------------------
agal_mafs_filt <- agal_mafs_filt %>% arrange(chromo, position)
cver_mafs_filt <- cver_mafs_filt %>% arrange(chromo, position)

# -----------------------------
# STEP 3: Sanity check
# -----------------------------
if(!all(agal_mafs_filt$position == cver_mafs_filt$position)) stop("Positions do not match!")
if(!all(agal_mafs_filt$chromo == cver_mafs_filt$chromo)) stop("Chromosomes do not match!")

cat("Shared SNPs:", nrow(agal_mafs_filt), "\n")
#Shared SNPs: 1,354,059 

# -----------------------------
# STEP 4: Define parameters
# -----------------------------
Ne <- c(1221, 1221)
gen <- c(0, 63)

# -----------------------------
# STEP 5: Build matrices
# -----------------------------
freq_mat <- data.frame(
  agal = agal_mafs_filt$knownEM,
  cver = cver_mafs_filt$knownEM
)

cov_mat <- data.frame(
  agal = agal_mafs_filt$nInd,
  cver = cver_mafs_filt$nInd
)

# -----------------------------
# STEP 6: Run chi-square test
# -----------------------------
pval <- adapted.chisq.test(
  freq = freq_mat,
  coverage = cov_mat,
  Ne = Ne,
  gen = gen
)

pval <- as.data.frame(pval)

# Check column name (important!)
print(colnames(pval))

# Adjust if needed:
pval$fdr <- p.adjust(pval[,1], method = "fdr")

# -----------------------------
# STEP 7: Attach metadata
# -----------------------------
freq <- cbind(pval, freq_mat)

freq$chromo <- agal_mafs_filt$chromo
freq$position <- agal_mafs_filt$position

# -----------------------------
# STEP 8: Extract SNPs under selection
# -----------------------------
freq_fdr <- subset(freq, fdr < 0.05)

freq_fdr$delta <- freq_fdr$cver - freq_fdr$agal

cat("SNPs under selection:", nrow(freq_fdr), "\n")
#SNPs under selection: 9,233 
cat("Percent selected:",
    round(nrow(freq_fdr)/nrow(freq)*100, 2), "%\n")
#ercent selected: 0.68 %
# -----------------------------
# STEP 9: Export selected SNPs
# -----------------------------
selected_snps <- freq_fdr %>%
  select(chromo, position)

write.table(
  selected_snps,
  file = "Psq_snp_list_sel.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

cat("Selection SNP list saved.\n")

# -----------------------------
# OPTIONAL: Neutral SNP extraction
# -----------------------------

# Create chrpos IDs
freq$chrpos <- paste(freq$chromo, freq$position, sep = "_")
freq_fdr$chrpos <- paste(freq_fdr$chromo, freq_fdr$position, sep = "_")

neutral_snps <- freq[!freq$chrpos %in% freq_fdr$chrpos, ]

cat("Neutral SNPs:", nrow(neutral_snps), "\n")
#Neutral SNPs: 1,344,826 

# Export neutral SNPs
write.table(
  neutral_snps[, c("chromo","position")],
  file = "Psq_snplist_neutral.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

cat("Neutral SNP list saved.\n")
