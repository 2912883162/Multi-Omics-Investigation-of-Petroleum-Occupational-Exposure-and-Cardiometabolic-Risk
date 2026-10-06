# ==============================================================================
# Script 04: Metagenomic Functional Pathway Stratification and Guild Contribution
# Description: Quantifies guild-specific contributions to endotoxin biosynthesis 
#              pathways (Kdo2-Lipid A and LPS superpathway) from stratified 
#              HUMAnN3 outputs.
# Figure: Figure 4A, Figure S8
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# ------------------------------------------------------------------------------
# 1. Target Pathway Definitions and HUMAnN3 Parsing
# ------------------------------------------------------------------------------

# Endotoxin biosynthesis pathways identified in the manuscript (Figure 5A-B)
TARGET_PATHWAYS <- list(
  "Lipid_A_Synthesis" = "KDO-NAGLIPASYN-PWY: superpathway of (Kdo)2-lipid A biosynthesis",
  "LPS_Biosynthesis"  = "LPSSYN-PWY: superpathway of lipopolysaccharide biosynthesis"
)

# Calculate relative guild contribution ratio to a specific functional pathway
calculate_guild_contribution <- function(path_abund_df, guild_membership, target_pwy_name, target_guild = "Guild_1") {
  guild_membership$Clean_Species <- str_replace(guild_membership$Species, "^s__", "")
  spec_map <- setNames(guild_membership$Guild_ID, guild_membership$Clean_Species)
  
  is_target <- grepl(target_pwy_name, path_abund_df[[1]], fixed = TRUE)
  sub_df <- path_abund_df[is_target, ]
  
  if (nrow(sub_df) == 0) {
    warning("Pathway not found in the input table: ", target_pwy_name)
    return(NULL)
  }
  
  strat_idx <- grep("\\|", sub_df[[1]])
  total_idx <- which(!grepl("\\|", sub_df[[1]]))
  sample_names <- colnames(path_abund_df)[-1]
  
  # Total community abundance for this pathway
  total_abundance <- as.numeric(sub_df[total_idx[1], sample_names])
  names(total_abundance) <- sample_names
  
  # Guild-specific stratified abundance
  strat_rows <- sub_df[strat_idx, ]
  species_extracted <- str_extract(strat_rows[[1]], "s__.*$") %>% str_replace("^s__", "")
  row_guilds <- spec_map[species_extracted]
  
  g_idx <- which(row_guilds == target_guild)
  if (length(g_idx) > 0) {
    guild_abundance <- colSums(as.matrix(strat_rows[g_idx, sample_names, drop = FALSE]))
  } else {
    guild_abundance <- rep(0, length(sample_names))
    names(guild_abundance) <- sample_names
  }
  
  # Contribution Ratio = Guild abundance / Total community pathway abundance
  res_df <- data.frame(
    Sample = sample_names,
    Pathway = target_pwy_name,
    Target_Guild = target_guild,
    Total_Abundance = total_abundance,
    Guild_Abundance = guild_abundance,
    Contribution_Ratio = ifelse(total_abundance > 1e-9, guild_abundance / total_abundance, 0)
  )
  
  return(res_df)
}

# ------------------------------------------------------------------------------
# 2. Between-Group Comparison of Guild Pathway Contributions
# ------------------------------------------------------------------------------

# Compare guild contribution ratios between occupational exposure groups
compare_guild_contributions <- function(contrib_df, metadata_df) {
  merged <- inner_join(contrib_df, metadata_df, by = "Sample")
  
  fit_model <- lm(Contribution_Ratio ~ Group, data = merged)
  coef_summary <- summary(fit_model)$coefficients
  
  message("Group difference test for: ", unique(contrib_df$Pathway))
  print(coef_summary)
  
  return(list(model = fit_model, summary = coef_summary))
}
