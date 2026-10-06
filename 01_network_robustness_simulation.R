# ==============================================================================
# Script 01: Ecological Network Robustness and Vulnerability Simulation
# Description: Evaluates topological robustness of microbial co-occurrence 
#              networks under random perturbations and degree-targeted attacks 
#              using Natural Connectivity (Yuan et al.).
# Figure: Figure 2F-G
# ==============================================================================

suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 1. Core Topological Metrics: Natural Connectivity and Network Attacks
# ------------------------------------------------------------------------------

# Calculate Natural Connectivity based on the adjacency matrix spectrum
calc_natural_connectivity <- function(g) {
  if (ecount(g) == 0 || vcount(g) <= 1) return(0)
  adj <- as_adjacency_matrix(g, sparse = FALSE)
  evals <- eigen(adj, only.values = TRUE)$values
  log(mean(exp(evals)))
}

# Simulate network degradation under random node removal
simulate_random_robustness <- function(g, removal_ratios = seq(0, 0.9, by = 0.05), n_sim = 50) {
  if (ecount(g) == 0) {
    return(data.frame(Removal_Ratio = removal_ratios, Connectivity = 0, SD = 0))
  }
  
  n_nodes <- vcount(g)
  res_mat <- matrix(0, nrow = n_sim, ncol = length(removal_ratios))
  
  for (i in seq_len(n_sim)) {
    shuffled_nodes <- sample(V(g)$name)
    for (j in seq_along(removal_ratios)) {
      n_rem <- floor(n_nodes * removal_ratios[j])
      if (n_rem == 0) {
        res_mat[i, j] <- calc_natural_connectivity(g)
      } else if (n_rem >= n_nodes) {
        res_mat[i, j] <- 0
      } else {
        g_sub <- delete_vertices(g, shuffled_nodes[1:n_rem])
        res_mat[i, j] <- calc_natural_connectivity(g_sub)
      }
    }
  }
  
  data.frame(
    Removal_Ratio = removal_ratios,
    Connectivity = colMeans(res_mat),
    SD = apply(res_mat, 2, sd),
    Simulation_Type = "Random Perturbation"
  )
}

# Simulate network degradation under degree-ranked hub removal
simulate_targeted_attack <- function(g, removal_ratios = seq(0, 0.5, by = 0.05)) {
  if (ecount(g) == 0) {
    return(data.frame(Removal_Ratio = removal_ratios, Connectivity = 0))
  }
  
  n_nodes <- vcount(g)
  deg <- degree(g)
  target_ranked_nodes <- names(sort(deg, decreasing = TRUE))
  
  conn_values <- numeric(length(removal_ratios))
  for (j in seq_along(removal_ratios)) {
    n_rem <- floor(n_nodes * removal_ratios[j])
    if (n_rem == 0) {
      conn_values[j] <- calc_natural_connectivity(g)
    } else if (n_rem >= n_nodes) {
      conn_values[j] <- 0
    } else {
      g_sub <- delete_vertices(g, target_ranked_nodes[1:n_rem])
      conn_values[j] <- calc_natural_connectivity(g_sub)
    }
  }
  
  data.frame(
    Removal_Ratio = removal_ratios,
    Connectivity = conn_values,
    Simulation_Type = "Targeted Attack"
  )
}

# Calculate Area Under the Curve (AUC) via trapezoidal numerical integration
calc_auc <- function(ratios, values) {
  sum(diff(ratios) * (head(values, -1) + tail(values, -1)) / 2)
}

# ------------------------------------------------------------------------------
# 2. Pipeline Execution: Guild Topological Robustness Profiling
# ------------------------------------------------------------------------------

run_guild_robustness_pipeline <- function(guild_graphs) {
  results_list <- list()
  auc_summary <- data.frame()
  
  for (guild_id in names(guild_graphs)) {
    g <- guild_graphs[[guild_id]]
    message("Evaluating network stability for: ", guild_id)
    
    # 1. Random perturbation simulation (0 - 90% node removal)
    res_random <- simulate_random_robustness(g, removal_ratios = seq(0, 0.9, by = 0.05), n_sim = 50)
    auc_rand <- calc_auc(res_random$Removal_Ratio, res_random$Connectivity)
    
    # 2. Targeted attack simulation (0 - 50% node removal)
    res_targeted <- simulate_targeted_attack(g, removal_ratios = seq(0, 0.5, by = 0.05))
    auc_targ <- calc_auc(res_targeted$Removal_Ratio, res_targeted$Connectivity)
    
    results_list[[guild_id]] <- list(random = res_random, targeted = res_targeted)
    auc_summary <- rbind(auc_summary, data.frame(
      Guild = guild_id,
      Random_AUC = round(auc_rand, 4),
      Targeted_AUC = round(auc_targ, 4)
    ))
  }
  
  return(list(curves = results_list, summary = auc_summary))
}
