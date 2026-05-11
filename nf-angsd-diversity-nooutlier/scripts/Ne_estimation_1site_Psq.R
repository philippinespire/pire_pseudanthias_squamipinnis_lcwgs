#####################################################
#######Ne estimation & Bootstrap for 1 Site#########
#####################################################
#Step 0: Load required libraries
setwd('/archive/carpenterlab/pire/pire_pseudanthias_squamipinnis_lcwgs/nf-angsd-diversity-nooutlier/results_all/angsd_pop/') #change based on your directory pathway

library(data.table)
library(boot)
library(tidyverse)
install.packages('R.utils')

# Step 1: Read in allele frequency data (*.maf.gz files)
#All SNPs
agal_mafs <- fread("PsqAGal_historic.mafs.gz", header = TRUE)
cver_mafs <- fread("PsqCVer_modern.mafs.gz", header = TRUE)

#neutral mafs
#agal_mafs <- fread("agal_sites_neutral.mafs.gz", header = TRUE)
#cver_mafs <- fread("cver_sites_neutral.mafs.gz", header = TRUE)

#Step 2: Rename columns for consistency
rename_mafs <- function(mafs) {
  setnames(mafs, c("knownEM", "nInd"), c("freq1", "nInd1"))
  return(mafs)
}

agal_mafs <- rename_mafs(agal_mafs)
cver_mafs <- rename_mafs(cver_mafs)

#Step 3: Merge populations for pairwise comparisons so that temporal replicates can be compared at each site
agal_cver_mafs <- agal_mafs[cver_mafs, on = .(chromo, position), 
                            .(chromo, position, 
                              freq1 = freq1, 
                              freq2 = i.freq1, 
                              nInd1 = nInd1, 
                              nInd2 = i.nInd1)]

#Step 4: Filter SNPs with maf > 0.07
agal_cver_mafs <- agal_cver_mafs[!is.na(freq1) & !is.na(freq2) & freq1 > 0.07 & freq2 > 0.07, ]


#Step 5: Define the Jorde & Ryman function for Ne estimation
jrNe2 <- function(maf1, maf2, n1, n2, gen) {
  Fsnum <- (maf1 - maf2)^2 + (1 - maf1 - (1 - maf2))^2
  z <- (maf1 + maf2) / 2
  z2 <- ((1 - maf1) + (1 - maf2)) / 2
  Fsdenom <- z * (1 - z) + z2 * (1 - z2)
  Fs <- sum(Fsnum, na.rm = TRUE) / sum(Fsdenom, na.rm = TRUE)
  
  sl <- 2 / (1 / n1 + 1 / n2)
  S <- length(maf1) * 2 / sum(2 / sl, na.rm = TRUE)
  S2 <- length(maf2) * 2 / sum(2 / n2, na.rm = TRUE)
  
  Fsprime <- (Fs * (1 - 1 / (4 * S)) - 1 / S) / ((1 + Fs / 4) * (1 - 1 / (2 * S2)))
  
  Ne <- gen / (2 * Fsprime)
  if (Ne < 0) Ne <- Inf  
  return(Ne)
}

#Step 6: Define the bootstrap function
jrNe2_boot <- function(data, indices, gen) {
  sampled_data <- data[indices, ]  # Resample rows
  return(jrNe2(sampled_data$freq1, sampled_data$freq2, sampled_data$nInd1, sampled_data$nInd2, gen))
}

#Step 7: Function to run Ne estimation and bootstrap confidence intervals per site
analyze <- function(data, gen, R = 1000) {
  Ne_estimate <- jrNe2(data$freq1, data$freq2, data$nInd1, data$nInd2, gen)
  
  boot_result <- boot(data = data, 
                      statistic = function(data, indices) jrNe2_boot(data, indices, gen), 
                      R = R)
  ci <- boot.ci(boot_result, type = "perc")
  
  CI_lower <- ifelse(is.null(ci$percent), NA, ci$percent[4])
  CI_upper <- ifelse(is.null(ci$percent), NA, ci$percent[5])
  
  return(list(Ne_estimate = Ne_estimate, CI_lower = CI_lower, CI_upper = CI_upper))
}

#Generation time of Pseudanthias squamipinnis is 1.8 years based on Fishbase, therefore, we divide the total number of years between historical & contemporary sampling periods by the generation time
#Albatross was collected in 1909 & 2022; 62.77 generations
#63 generations was used
result_agal_cver <- analyze(agal_cver_mafs, 63)

#Step 8: Print results
cat(sprintf("Ne for AGal vs CVer (Site 1): %.2f (95%% CI: %.2f - %.2f)\n", 
            result_agal_cver$Ne_estimate, result_agal_cver$CI_lower, result_agal_cver$CI_upper))

#############################
############RESULTS##########
#############################
##All SNPs Ne for AGal vs CVer (Site 1): 1450.34 (95% CI: 1423.68 - 1479.49) for 0.08 MAF
##All SNPs Ne or AGal vs CVer (Site 1): 1422.56 (95% CI: 1396.49 - 1446.93) for 0.07 MAF
##Neutral 