############Directory Set-up##################
setwd("/archive/carpenterlab/pire/pire_pseudanthias_squamipinnis_lcwgs/nf-angsd-diversity-nooutlier/plots")

# ALL SNPs
dir_all_theta  <- "/archive/carpenterlab/pire/pire_pseudanthias_squamipinnis_lcwgs/nf-angsd-diversity-nooutlier/results_all/angsd_pop_theta"
dir_all_depth  <- "/archive/carpenterlab/pire/pire_pseudanthias_squamipinnis_lcwgs/nf-angsd-diversity-nooutlier/results_all/depth"

# NEUTRAL SNPs
dir_neut_theta <- "/archive/carpenterlab/pire/pire_pseudanthias_squamipinnis_lcwgs/nf-angsd-diversity-nooutlier/results/angsd_pop_theta"
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
