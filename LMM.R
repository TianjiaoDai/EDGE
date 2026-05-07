# Linear Mixed-effects Models for Diversity Analysis
# Performs LMM to test Treatment effects within each site, accounting for Block as a random factor in the EDGE network.

library(lme4)
library(lmerTest) 
library(car)
library(tidyverse)

# 1. Define core analysis function 
run_site_lmm <- function(div_data, meta_data, site_name) {
  message(paste("Processing Site:", site_name, "at", Sys.time()))
  
  # Loop through each diversity index
  results_list <- lapply(colnames(div_data), function(idx_name) {
    
    # Integrate data frames
    df_temp <- data.frame(
      Value = div_data[[idx_name]],
      Treatment = factor(meta_data$Treatment),
      Block = factor(meta_data$Block)
    )
    
    # Fit linear mixed-effects model
    # Random effect (1|Block) accounts for the randomized block design
    model <- lmer(Value ~ Treatment + (1|Block), data = df_temp)
    
    # Extract statistics (Type II Wald chi-square tests)
    anova_res <- car::Anova(model, type = 2)
    sum_mod <- summary(model)
    coef_table <- coef(sum_mod)
    
    # Extract key metrics: Estimate, SE, t-value, and P-value
    res <- c(
      Estimate  = coef_table["TreatmentDrought", "Estimate"],
      Std.Error = coef_table["TreatmentDrought", "Std. Error"],
      t.value   = coef_table["TreatmentDrought", "t value"],
      Chisq     = anova_res["Treatment", "Chisq"],
      P.value   = anova_res["Treatment", "Pr(>Chisq)"]
    )
    return(res)
  })
  
  # Organize output into a data frame
  res_matrix <- do.call(rbind, results_list)
  rownames(res_matrix) <- colnames(div_data)
  return(as.data.frame(res_matrix))
}

# 2. Batch process all sites (Continental-scale analysis) 
# divindex: diversity index matrix; treatused: metadata
site_list <- unique(treatused$Site)

all_site_results <- lapply(site_list, function(site) {
  # Subset data for the specific site
  site_idx  <- treatused$Site == site
  curr_div  <- divindex[site_idx, ]
  curr_meta <- treatused[site_idx, ]
  
  # Run LMM for the current site
  site_res <- run_site_lmm(curr_div, curr_meta, site)
  
  # Save results for each site
  write.csv(site_res, paste0(site, "_16S_div_compare_Chronic.csv"))
  
  return(site_res)
})

names(all_site_results) <- site_list