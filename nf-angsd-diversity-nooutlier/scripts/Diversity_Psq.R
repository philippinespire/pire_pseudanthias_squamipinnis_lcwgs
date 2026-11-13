############Directory Set-up##################
setwd("/archive/carpenterlab/pire/pire_pseudanthias_squamipinnis_lcwgs/nf-angsd-diversity-nooutlier/plots")

# ALL SNPs
dir_all_theta  <- "/archive/carpenterlab/pire/pire_pseudanthias_squamipinnis_lcwgs/nf-angsd-diversity-nooutlier/results_all/angsd_pop_theta"
dir_all_depth  <- "/archive/carpenterlab/pire/pire_pseudanthias_squamipinnis_lcwgs/nf-angsd-diversity-nooutlier/results_all/depth"

# NEUTRAL SNPs
dir_neut_theta <- "/archive/carpenterlab/pire/pire_pseudanthias_squamipinnis_lcwgs/nf-angsd-diversity-nooutlier/results_neutral/angsd_pop_theta"
dir_neut_depth <- "/archive/carpenterlab/pire/pire_pseudanthias_squamipinnis_lcwgs/nf-angsd-diversity-nooutlier/results_all/depth"

############Load Packages####################
library(tidyverse)
library(cowplot)
library(boot)
library(purrr)

#############Functions#################
load_theta <- function(dir, file, era, dataset){
  read_table(file.path(dir, file)) %>%
    mutate(
      Era = era,
      Dataset = dataset
    )
}

process_depth <- function(dir, file, n_ind){
  read_table(file.path(dir, file)) %>%
    group_by(chr) %>%
    summarise(avg_Depth = mean(totDepth), .groups="drop") %>%
    mutate(ind_depth = avg_Depth / n_ind) %>%
    rename(Chr = chr)
}

############Filter based on depth#################
filter_depth_auto <- function(df, method="quantile", lower=0.05, upper=0.95){
  
  df <- df %>% filter(is.finite(ind_depth))
  
  if(method == "mean_sd"){
    mu <- mean(df$ind_depth, na.rm=TRUE)
    sdv <- sd(df$ind_depth, na.rm=TRUE)
    min_val <- mu - sdv
    max_val <- mu + sdv
    
  } else if(method == "quantile"){
    min_val <- quantile(df$ind_depth, lower, na.rm=TRUE)
    max_val <- quantile(df$ind_depth, upper, na.rm=TRUE)
    
  } else {
    stop("Method must be 'mean_sd' or 'quantile'")
  }
  
  cat("Depth filter range:", round(min_val,3), "-", round(max_val,3), "\n")
  
  df %>% filter(ind_depth >= min_val & ind_depth <= max_val)
}

############Load Files#####################
# ==== All SNPs ====
hist_all <- load_theta(dir_all_theta, "PsqAGal_historic.pestPG", "Historical", "All SNPs")
mod_all  <- load_theta(dir_all_theta, "PsqCVer_modern.pestPG", "Contemporary", "All SNPs")

# ==== Neutral SNPs ====
hist_neut <- load_theta(dir_neut_theta, "PsqAGal_historic.pestPG", "Historical", "Neutral SNPs")
mod_neut  <- load_theta(dir_neut_theta, "PsqCVer_modern.pestPG", "Contemporary", "Neutral SNPs")

# Combine
theta_raw <- bind_rows(hist_all, mod_all, hist_neut, mod_neut)

############Depth#################
depth <- process_depth(
  dir_all_depth,
  "depth_all_sites.pos.gz",
  n_ind = 62
)

############Merge & metrics########
theta <- theta_raw %>%
  left_join(depth, by="Chr") %>%
  mutate(
    tW_bysite = tW / nSites,
    tP_bysite = tP / nSites,
    Era = factor(Era, levels=c("Historical","Contemporary")),
    Dataset = factor(Dataset, levels=c("All SNPs","Neutral SNPs"))
  )

############ APPLY FILTER ###########
theta_filt <- filter_depth_auto(theta, method="quantile", lower=0.05, upper=0.95)
#Depth filter range: Depth filter range: 0.705 - 1.161 

############Bootstrap#############
boot_metric <- function(x, R=1000){
  x <- na.omit(x)
  if(length(x) < 2){
    return(tibble(mean=NA, se=NA))
  }
  b <- boot(x, function(x,i) mean(x[i]), R=R)
  tibble(
    mean = mean(x),
    se   = sd(b$t)
  )
}

boot_results <- theta_filt %>%
  group_by(Dataset, Era) %>%
  summarise(
    pi = list(boot_metric(tP_bysite)),
    theta = list(boot_metric(tW_bysite)),
    .groups="drop"
  ) %>%
  unnest(cols=c(pi, theta), names_sep="_")

############ PREP DATA ############
plot_long <- boot_results %>%
  pivot_longer(
    cols = c(pi_mean, theta_mean, pi_se, theta_se),
    names_to = c("Metric", ".value"),
    names_pattern = "(pi|theta)_(mean|se)"
  ) %>%
  mutate(
    Metric = recode(Metric,
                    pi = "Pi (π)",
                    theta = "Theta (θ)")
  )

############ FINAL PLOT ############
p_final <- ggplot(plot_long,
                  aes(x=Era,
                      y=mean,
                      group=interaction(Dataset, Metric),
                      color=Dataset,
                      shape=Metric)) +
  
  geom_line(size=1, alpha=0.8,
            position=position_dodge(width=0.15)) +
  
  geom_point(size=3,
             position=position_dodge(width=0.15)) +
  
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=0.08,
                alpha=0.7,
                position=position_dodge(width=0.15)) +
  
  facet_grid(Metric ~ Dataset, scales="free_y") +
  
  theme_classic() +
  
  labs(
    title = "Genetic Diversity Change Over Time",
    subtitle = "Pseudanthias squamipinnis",
    x = "",
    y = "Diversity (per site)",
    color = "Dataset",
    shape = "Metric"
  ) +
  
  scale_shape_manual(values=c(16, 17)) +
  
  theme(
    strip.text = element_text(face="bold"),
    plot.title = element_text(face="bold"),
    panel.spacing = unit(1, "lines")
  )

print(p_final)

ggsave(
  "pi_theta_plot.pdf",   # vector format (best)
  p_final,
  width = 10,
  height = 8
)

ggsave(
  "pi_theta_plot.png",   # high-res raster backup
  p_final,
  width = 10,
  height = 8,
  dpi = 600
)

############ ########################### ############
############ Test for Significance: PI & THETA ######
############ ########################### ############

test_metric_shift <- function(df, metric){
  
  df_wide <- df %>%
    select(Chr, Dataset, Era, all_of(metric)) %>%
    pivot_wider(names_from = Era, values_from = all_of(metric)) %>%
    mutate(delta = Contemporary - Historical)
  
  df_wide %>%
    group_by(Dataset) %>%
    summarise(
      mean_delta = mean(delta, na.rm=TRUE),
      t_p = t.test(delta, mu=0)$p.value,
      wilcox_p = wilcox.test(delta, mu=0)$p.value,
      .groups="drop"
    )
}

pi_tests <- test_metric_shift(theta_filt, "tP_bysite")
theta_tests <- test_metric_shift(theta_filt, "tW_bysite")

print(pi_tests)
#Dataset      mean_delta   t_p  wilcox_p
# 1 All SNPs         0.00344 2.07e- 66 1.27e-90
#2 Neutral SNPs   -0.00921 2.71e-313 0

print(theta_tests)
#Dataset      mean_delta   t_p wilcox_p
#  1 All SNPs         -0.0354     0        0
#2 Neutral SNPs    -0.0437     0        0
################################################
################### TAJIMA'S D ##################
################################################

plot_tajima <- function(df){
  ggplot(df, aes(x=Tajima, fill=Era)) +
    geom_density(alpha=0.6) +
    facet_wrap(~Dataset) +   
    theme_classic() +
    labs(title="Tajima's D Distribution", x="Tajima's D", y="Density")
}

print(plot_tajima(theta_filt))

######## Tajima's D analysis ###################

tajima_shift <- theta_filt %>%
  select(Chr, Dataset, Era, Tajima) %>%
  pivot_wider(names_from = Era, values_from = Tajima) %>%
  mutate(delta_Tajima = Contemporary - Historical)

tajima_summary <- tajima_shift %>%
  group_by(Dataset) %>%
  summarise(
    mean_delta = mean(delta_Tajima, na.rm=TRUE),
    median_delta = median(delta_Tajima, na.rm=TRUE),
    prop_positive = mean(delta_Tajima > 0, na.rm=TRUE),
    prop_negative = mean(delta_Tajima < 0, na.rm=TRUE),
    .groups="drop"
  )

tajima_tests <- tajima_shift %>%
  group_by(Dataset) %>%
  summarise(
    p_value = wilcox.test(delta_Tajima, mu = 0)$p.value,
    .groups="drop"
  )

print(tajima_summary)
#Dataset      mean_delta median_delta prop_positive prop_negative
#  1 All SNPs          0.665        0.654         0.992       0.00842
#2 Neutral SNPs      0.617        0.620         0.973       0.0271 

print(tajima_tests)
#Dataset        p_value
#  1 All SNPs     0
#2 Neutral SNPs 0

################################################
######## TEMPORAL DIFFERENTIATION ###############
################################################

temporal_diff <- tajima_shift %>%
  group_by(Dataset) %>%
  summarise(
    mean_abs_change = mean(abs(delta_Tajima), na.rm=TRUE),
    prop_large_shift = mean(abs(delta_Tajima) > 0.5, na.rm=TRUE),
    .groups="drop"
  )

print(temporal_diff)
#Dataset      mean_abs_change prop_large_shift
#  1 All SNPs               0.668            0.753
#2 Neutral SNPs           0.625            0.667

################################################
######## ALLELE FREQUENCY CHANGE (Δp) ##########
################################################

delta_freq <- theta_filt %>%
  select(Chr, Dataset, Era, tP_bysite) %>%
  pivot_wider(names_from = Era, values_from = tP_bysite) %>%
  mutate(delta = Contemporary - Historical)

ggplot(delta_freq, aes(x=delta, fill=Dataset)) +
  geom_histogram(alpha=0.6, bins=50, position="identity") +
  facet_wrap(~Dataset) +
  theme_classic() +
  labs(
    title="Distribution of Allele Frequency Change (Δp)",
    x="Change (Contemporary - Historical)",
    y="Count"
  )

################################################
############ TEMPORAL FST ANALYSIS #############
################################################

############ LOAD LIBRARIES ############
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)

############ SET DIRECTORIES ############
dir_all <- "/archive/carpenterlab/pire/pire_pseudanthias_squamipinnis_lcwgs/nf-angsd-diversity-nooutlier/results_all/angsd_pop/"
dir_neutral <- "/archive/carpenterlab/pire/pire_pseudanthias_squamipinnis_lcwgs/nf-angsd-diversity-nooutlier/results_neutral/angsd_pop/"

############ FUNCTION TO READ FST ############
read_fst_windowed <- function(file, dataset_label){
  
  lines <- readLines(file)
  lines <- lines[lines != ""]
  
  # Remove header line
  lines <- lines[!grepl("^region\\s+chr", lines)]
  
  # Split by whitespace
  split_lines <- strsplit(lines, "\\s+")
  max_cols <- max(lengths(split_lines))
  
  split_padded <- lapply(split_lines, function(x){
    length(x) <- max_cols
    return(x)
  })
  
  df <- as.data.frame(do.call(rbind, split_padded), stringsAsFactors = FALSE)
  
  # Clean bad values
  df[df %in% c("nan", "inf", "-inf")] <- NA
  
  # Convert types
  df[] <- lapply(df, type.convert, as.is = TRUE)
  
  # Assign columns
  if(ncol(df) == 4){
    
    colnames(df) <- c("chr","midPos","Nsites","FST")
    
    df <- df %>%
      mutate(midPos = as.numeric(midPos))
    
  } else if(ncol(df) == 5){
    
    colnames(df) <- c("chr","start","end","Nsites","FST")
    
    df <- df %>%
      mutate(
        start  = as.numeric(start),
        end    = as.numeric(end),
        midPos = (start + end) / 2
      )
    
  } else {
    stop(paste("Unexpected format in", file))
  }
  
  # Final cleanup
  df %>%
    mutate(
      FST  = as.numeric(FST),
      Pair = gsub(".windowed.fst.txt", "", basename(file)),
      Dataset = dataset_label
    ) %>%
    filter(is.finite(FST))
}

############ LOAD FILES (BOTH DATASETS) ############
fst_files_all <- list.files(dir_all, pattern = "windowed.fst.txt$", full.names = TRUE)
fst_files_neutral <- list.files(dir_neutral, pattern = "windowed.fst.txt$", full.names = TRUE)

############ READ ALL DATA ############
fst_data_all <- map_dfr(fst_files_all, ~read_fst_windowed(.x, "All"))
fst_data_neutral <- map_dfr(fst_files_neutral, ~read_fst_windowed(.x, "Neutral"))

############ COMBINE ############
fst_data <- bind_rows(fst_data_all, fst_data_neutral)

############ CHECK OUTPUT ############
head(fst_data)
summary(fst_data$FST)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-0.001804  0.013250  0.016018  0.016914  0.019469  0.078021 

################################################
############ GLOBAL FST SUMMARY #################
################################################

fst_summary <- fst_data %>%
  group_by(Dataset, Pair) %>%
  summarise(
    mean_fst = mean(FST, na.rm=TRUE),
    median_fst = median(FST, na.rm=TRUE),
    lower_ci = quantile(FST, 0.025, na.rm=TRUE),
    upper_ci = quantile(FST, 0.975, na.rm=TRUE),
    p_wilcox = wilcox.test(FST, mu = 0)$p.value,
    prop_gt0 = mean(FST > 0, na.rm=TRUE),
    .groups="drop"
  )

print(fst_summary)

# Dataset Pair                               mean_fst median_fst lower_ci upper_ci p_wilcox prop_gt0
#  1 All     PsqAGal_historic_vs_PbbCVer_modern   0.0165     0.0155  0.00836   0.0310        0     1.00
#2 Neutral PsqAGal_historic_vs_PsqCVer_modern   0.0173     0.0165  0.00905   0.0301        0     1.00

################################################
############ MANHATTAN-STYLE PLOT ##############
################################################

fst_data_clean <- fst_data %>%
  filter(!is.na(FST), is.finite(FST))

ggplot(fst_data_clean, aes(x=midPos, y=FST, color=Dataset)) +
  geom_point(alpha=0.5, size=0.6) +
  facet_wrap(~Pair, scales="free_x") +
  theme_bw() +
  labs(
    title="Temporal FST Across Genome (All vs Neutral SNPs)",
    x="Genomic Position",
    y="FST"
  )

################################################
############ OUTLIER WINDOWS ###################
################################################

fst_outliers <- fst_data_clean %>%
  group_by(Dataset, Pair) %>%
  mutate(threshold = quantile(FST, 0.99, na.rm=TRUE)) %>%
  filter(FST >= threshold)

cat("Number of high-FST windows:\n")
print(fst_outliers %>% count(Dataset, Pair))
#Dataset Pair                                   n
#  1 All     PbbAGal_historic_vs_PbbCVer_modern    79
#2 Neutral PbbAGal_historic_vs_PbbCVer_modern    79

################################################
############### SAVE FST OUTPUTS ###############
################################################
write.csv(fst_summary, "fst_summary_all_vs_neutral.csv", row.names=FALSE)
write.csv(fst_outliers, "fst_outlier_windows_all_vs_neutral.csv", row.names=FALSE)
