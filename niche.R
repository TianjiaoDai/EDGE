# Niche Analysis
# Comprehensive workflow: pCCA niche axis extraction, taxon-level niche quantification, distribution testing, and community-weighted mean (CWM) calculation.

library(vegan)
library(tidyverse)
library(stats)

# 1. pCCA: Extracting the Niche Axis -----------------------------------------
# Synthetic environmental axes are generated via PCA (excluding moisture)
# Moisture (SWC) is then conditioned out to isolate orthogonal niche axes

pca_env <- prcomp(envs[sampid, -1], center = TRUE, scale. = TRUE)
env_pca_scores <- as.data.frame(pca_env$x[, 1:2]) 

# partial CCA to obtain environmental coordinates independent of SWC
cca_data <- data.frame(water = envs[sampid, 1], env_pca_scores)
pcca_mod <- cca(commi ~ PC1 + PC2 + Condition(water), data = cca_data)

# Extract sample coordinates along the primary niche axis (CCA1)
niche_axis <- scores(pcca_mod, display = "lc", choices = 1)

# 2. Taxon-level Niche Quantification (Eq 2 & 3) -----------------------------
# Function to calculate niche optimum and width for each taxon

get_taxon_niche <- function(comm, axis_coord) {
  # Calculate relative abundance (b/sum(b))
  rel_abund <- apply(comm, 2, function(x) x / sum(x))
  
  # Eq 2: Niche Optimum (abundance-weighted mean)
  optimum <- apply(rel_abund, 2, function(b) sum(b * axis_coord))
  
  # Eq 3: Niche Width (abundance-weighted standard deviation)
  # Standardized weighted SD logic
  width <- apply(comm, 2, function(b) {
    w_mean <- weighted.mean(axis_coord, b)
    sqrt(sum(b * (axis_coord - w_mean)^2) / sum(b))
  })
  
  return(data.frame(Taxon = colnames(comm), Optimum = optimum, Width = width))
}

# Quantify niches for Control and Chronic treatments separately
niche_control <- get_taxon_niche(commi[treat$Treatment == "Control", ], 
                                niche_axis[treat$Treatment == "Control"])
niche_chronic <- get_taxon_niche(commi[treat$Treatment == "Chronic", ], 
                                niche_axis[treat$Treatment == "Chronic"])

# 3. Statistical Comparison (K-W Test) ---------------------------------------
# Comparing the probability distributions of niche attributes between groups

niche_combined <- bind_rows(
  niche_control %>% mutate(Group = "Control"),
  niche_chronic %>% mutate(Group = "Chronic")
)

kw_optimum <- kruskal.test(Optimum ~ factor(Group), data = niche_combined)
kw_width   <- kruskal.test(Width ~ factor(Group), data = niche_combined)

message(paste("K-W Test (Optimum) P-value:", round(kw_optimum$p.value, 5)))
message(paste("K-W Test (Width) P-value:", round(kw_width$p.value, 5)))

# 4. Community-weighted Mean (CWM) Niche -------------------------------------
# Scaling taxon-level attributes to the whole community level

calc_cwm <- function(comm_matrix, taxa_niche_attr) {
  # Standardize relative abundance per sample
  rel_comm <- comm_matrix / rowSums(comm_matrix)
  
  # Alignment check for taxa names
  common_taxa <- intersect(colnames(rel_comm), names(taxa_niche_attr))
  rel_comm <- rel_comm[, common_taxa]
  niche_vec <- taxa_niche_attr[common_taxa]
  
  # Matrix multiplication for weighted average
  cwm_values <- as.matrix(rel_comm) %*% as.matrix(niche_vec)
  return(as.data.frame(cwm_values))
}

# Calculate CWM Niche Width for each sample
# (Uses taxon widths calculated in Step 2)
taxon_widths <- setNames(niche_combined$Width, niche_combined$Taxon) # Simplified example
cwm_niche_width <- calc_cwm(commi, taxon_widths)