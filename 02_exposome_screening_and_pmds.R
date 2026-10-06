# ==============================================================================
# Script 02: Exposome Feature Screening and PMDS/GDI Index Construction
# Description: Two-stage feature selection (non-negative Elastic Net and Bayesian 
#              Kernel Machine Regression) on targeted plasma exposome, and 
#              construction of the Gut Dysbiosis Index (GDI) and Petroleum 
#              Mixture-Derived Signature (PMDS).
# Figure: Figure 3C-E, 3G
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(glmnet)
  library(bkmr)
})

# ------------------------------------------------------------------------------
# 1. Gut Dysbiosis Index (GDI) Calculation
# ------------------------------------------------------------------------------

# Composite scoring based on prioritized ecological guilds: (Z_G1 + Z_G3) - Z_G6
calculate_gdi_score <- function(guild_data) {
  z_g1 <- as.numeric(scale(guild_data$Guild_1))
  z_g3 <- as.numeric(scale(guild_data$Guild_3))
  z_g6 <- as.numeric(scale(guild_data$Guild_6))
  
  gdi <- (z_g1 + z_g3) - z_g6
  return(gdi)
}

# Rank-based Inverse Normal Transformation for exposure concentrations
quantile_normalisation <- function(x) {
  valid <- !is.na(x)
  qnorm((rank(x[valid], ties.method = "average") - 0.5) / sum(valid))
}

# ------------------------------------------------------------------------------
# 2. Stage 1: Elastic Net Feature Selection with Non-Negative Constraints
# ------------------------------------------------------------------------------

# Penalized regression to prioritize candidate exposures positively linked to GDI
# Specification: alpha = 0.5, lower.limits = 0, 10-fold cross-validation
screen_exposures_elastic_net <- function(exposure_matrix, gdi_vector, alpha = 0.5) {
  set.seed(123)
  cv_model <- cv.glmnet(
    x = as.matrix(exposure_matrix),
    y = as.numeric(gdi_vector),
    alpha = alpha,
    lower.limits = 0,
    nfolds = 10
  )
  
  # Extract non-zero positive coefficients at lambda.min
  model_coefs <- coef(cv_model, s = "lambda.min")
  retained <- data.frame(
    Variable = rownames(model_coefs),
    Coefficient = as.numeric(model_coefs)
  ) %>%
    filter(Variable != "(Intercept)", Coefficient > 0) %>%
    arrange(desc(Coefficient))
  
  return(list(fit = cv_model, retained_features = retained))
}

# ------------------------------------------------------------------------------
# 3. Stage 2: Bayesian Kernel Machine Regression (BKMR) Mixture Modeling
# ------------------------------------------------------------------------------

# Estimate joint mixture effects and Posterior Inclusion Probabilities (PIPs)
# Specification: 50,000 MCMC iterations, varsel = TRUE, PIP >= 0.5 cutoff
fit_bkmr_mixture_selection <- function(gdi_outcome, exposure_mixture, n_iter = 50000) {
  set.seed(1234)
  
  fit_km <- kmbayes(
    y = as.numeric(gdi_outcome),
    Z = as.matrix(exposure_mixture),
    iter = n_iter,
    verbose = FALSE,
    varsel = TRUE
  )
  
  pips <- ExtractPIPs(fit_km)
  retained_core <- pips %>% filter(PIP >= 0.5)
  
  return(list(bkmr_fit = fit_km, pips = pips, core_chemicals = retained_core))
}

# ------------------------------------------------------------------------------
# 4. Petroleum Mixture-Derived Signature (PMDS) Weighted Synthesis
# ------------------------------------------------------------------------------

# Construct composite PMDS based on verified relative component contributions
# Relative weights: 3-PBA (24.8%), PFOS (22.2%), 8:2 Cl-PFESA (21.3%), PNP (18.7%), DBP (12.9%)
calculate_pmds_score <- function(chem_data) {
  weights <- c(
    "3-PBA"        = 0.248,
    "PFOS"         = 0.222,
    "8:2 Cl-PFESA" = 0.213,
    "PNP"          = 0.187,
    "DBP"          = 0.129
  )
  
  z_chem <- scale(chem_data[, names(weights)])
  pmds <- as.numeric(z_chem %*% weights)
  return(pmds)
}
