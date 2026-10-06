# ==============================================================================
# Script 05: Proteomic WGCNA Module Detection and PPI Subnetwork Analysis
# Description: Implements standard Weighted Gene Co-expression Network Analysis 
#              (WGCNA) to identify co-expressed protein modules correlated with 
#              Guild 1 abundance and LDL-C, followed by extraction and topological 
#              characterization of the module-specific PPI subnetwork.
# Figure: Figure 4D-E
# ==============================================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(igraph)
  library(dplyr)
})

options(stringsAsFactors = FALSE)
try(enableWGCNAThreads(), silent = TRUE)

# ------------------------------------------------------------------------------
# 1. Standard WGCNA Network Construction and Module Identification
# ------------------------------------------------------------------------------

# Execute standard WGCNA pipeline
# Specification: Soft-thresholding R^2 >= 0.8, minModuleSize = 15, mergeCutHeight = 0.20
run_standard_wgcna <- function(protein_matrix, min_mod_size = 15, merge_cut = 0.20, r2_target = 0.80) {
  
  # Step 1: Select soft-thresholding power
  sft <- pickSoftThreshold(
    protein_matrix,
    powerVector = 1:20,
    networkType = "signed",
    verbose = 0
  )
  
  fit_df <- sft$fitIndices
  qualified <- fit_df$Power[!is.na(fit_df$SFT.R.sq) & fit_df$SFT.R.sq >= r2_target]
  selected_power <- if (length(qualified) > 0) min(qualified) else fit_df$Power[which.max(fit_df$SFT.R.sq)]
  
  message("Selected soft-thresholding power: ", selected_power)
  
  # Step 2: One-step blockwise network construction and module detection
  net <- blockwiseModules(
    protein_matrix,
    power = selected_power,
    networkType = "signed",
    TOMType = "signed",
    minModuleSize = min_mod_size,
    mergeCutHeight = merge_cut,
    numericLabels = FALSE,
    pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    verbose = 0
  )
  
  # Step 3: Extract Module Eigengenes (MEs)
  module_colors <- net$colors
  module_mes <- orderMEs(net$MEs)
  
  message("Identified co-expression modules: ", length(unique(module_colors)))
  
  return(list(
    power = selected_power,
    colors = module_colors,
    mes = module_mes,
    net_object = net
  ))
}

# ------------------------------------------------------------------------------
# 2. Module-Trait Relationship Correlation
# ------------------------------------------------------------------------------

# Correlate Module Eigengenes with Guild 1 abundance and LDL-C traits
correlate_modules_with_traits <- function(mes_df, traits_df) {
  me_cols <- colnames(mes_df)
  trait_cols <- colnames(traits_df)
  
  cor_mat <- matrix(NA, nrow = length(me_cols), ncol = length(trait_cols),
                    dimnames = list(me_cols, trait_cols))
  p_mat <- matrix(NA, nrow = length(me_cols), ncol = length(trait_cols),
                  dimnames = list(me_cols, trait_cols))
  
  for (i in seq_along(me_cols)) {
    for (j in seq_along(trait_cols)) {
      x <- as.numeric(mes_df[[me_cols[i]]])
      y <- as.numeric(traits_df[[trait_cols[j]]])
      keep <- complete.cases(x, y)
      
      if (sum(keep) >= 5) {
        test_res <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
        cor_mat[i, j] <- test_res$estimate
        p_mat[i, j] <- test_res$p.value
      }
    }
  }
  
  return(list(correlation = cor_mat, p_values = p_mat))
}

# ------------------------------------------------------------------------------
# 3. Objective Module PPI Subnetwork Extraction and Centrality Analysis
# ------------------------------------------------------------------------------

# Objectively extract the largest interacting subnetwork and rank topological hubs
# Specification: STRING medium confidence cutoff (combined_score >= 400)
build_module_ppi_network <- function(ppi_edge_list, module_proteins, min_score = 400) {
  
  # Filter PPI interactions where both nodes reside within the prioritized module
  filtered_edges <- ppi_edge_list %>%
    filter(
      from %in% module_proteins,
      to %in% module_proteins,
      combined_score >= min_score
    )
  
  if (nrow(filtered_edges) == 0) {
    warning("No interacting protein pairs found within the module at cutoff: ", min_score)
    return(NULL)
  }
  
  g <- graph_from_data_frame(filtered_edges, directed = FALSE)
  
  # Objectively extract the largest connected component
  components_info <- igraph::components(g)
  largest_comp_id <- which.max(components_info$csize)
  main_nodes <- names(components_info$membership[components_info$membership == largest_comp_id])
  g_sub <- induced_subgraph(g, main_nodes)
  
  # Calculate topological network centralities
  deg_vec <- degree(g_sub)
  bet_vec <- betweenness(g_sub)
  
  # Rank hub proteins by node degree
  hub_ranking <- data.frame(
    Protein = names(deg_vec),
    Degree = as.numeric(deg_vec),
    Betweenness = round(as.numeric(bet_vec), 2)
  ) %>%
    arrange(desc(Degree))
  
  message("Extracted core PPI subnetwork: Nodes = ", vcount(g_sub), ", Edges = ", ecount(g_sub))
  
  return(list(graph = g_sub, hub_table = hub_ranking))
}
