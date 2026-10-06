# ==============================================================================
# Script 03: Mediation Analysis and Restricted Cubic Spline (RCS) Modeling
# Description: Evaluates the intermediate roles of gut dysbiosis (GDI and Guild 1) 
#              connecting PMDS exposure to LDL-C elevation via causal mediation 
#              analysis, and fits dose-response curves via RCS.
# Figure: Figure 3I, Supplementary Figure S6
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(mediation)
  library(rms)
})

# ------------------------------------------------------------------------------
# 1. Restricted Cubic Spline (RCS) Modeling
# ------------------------------------------------------------------------------

# Fit dose-response curves for PMDS and LDL-C using restricted cubic splines
fit_rcs_exposure_response <- function(df, n_knots = 4) {
  old_dist <- options(datadist = "dd")
  .GlobalEnv$dd <- datadist(df)
  on.exit({
    options(old_dist)
    if (exists("dd", envir = .GlobalEnv)) rm("dd", envir = .GlobalEnv)
  })
  
  rcs_fit <- ols(LDL ~ rcs(PMDS, n_knots), data = df)
  rcs_anova <- anova(rcs_fit)
  
  return(list(model = rcs_fit, anova = rcs_anova))
}

# ------------------------------------------------------------------------------
# 2. Nonparametric Bootstrap Causal Mediation Analysis
# ------------------------------------------------------------------------------

# Test the mediation effect: PMDS (Exposure) -> Mediator (GDI or Guild 1) -> LDL (Outcome)
# Covariates can optionally be specified (e.g., covariates = c("Age", "Sex", "BMI"))
run_bootstrap_mediation <- function(analysis_data, mediator_name = "GDI_Score", covariates = NULL, n_sims = 1000) {
  set.seed(42)
  
  cov_str <- if (!is.null(covariates) && length(covariates) > 0) paste("+", paste(covariates, collapse = " + ")) else ""
  
  # Mediator model: PMDS -> Mediator
  form_med <- as.formula(paste(mediator_name, "~ PMDS", cov_str))
  fit_m <- eval(bquote(lm(.(form_med), data = analysis_data)))
  
  # Outcome model: PMDS + Mediator -> LDL
  form_out <- as.formula(paste("LDL ~ PMDS +", mediator_name, cov_str))
  fit_y <- eval(bquote(lm(.(form_out), data = analysis_data)))
  
  # Nonparametric bootstrap mediation (1,000 resamples)
  med_out <- mediate(
    model.m = fit_m,
    model.y = fit_y,
    treat = "PMDS",
    mediator = mediator_name,
    boot = TRUE,
    sims = n_sims
  )
  
  s <- summary(med_out)
  
  res_df <- data.frame(
    Mediator         = mediator_name,
    ACME_Estimate    = round(s$d0, 4),           # Average causal mediation effect
    ACME_95CI_Low    = round(s$d0.ci[1], 4),
    ACME_95CI_High   = round(s$d0.ci[2], 4),
    ACME_Pvalue      = round(s$d0.p, 4),
    ADE_Estimate     = round(s$z0, 4),           # Average direct effect
    ADE_Pvalue       = round(s$z0.p, 4),
    Total_Effect     = round(s$tau.coef, 4),     # Total effect
    Total_Pvalue     = round(s$tau.p, 4),
    Prop_Mediated    = round(s$n0 * 100, 2),     # Proportion mediated (%)
    Prop_Mediated_P  = round(s$n0.p, 4)
  )
  
  return(list(mediation_model = med_out, summary_table = res_df))
}

# ------------------------------------------------------------------------------
# 3. Pipeline Execution: Evaluating GDI and Guild 1 Pathways
# ------------------------------------------------------------------------------

run_full_mediation_pipeline <- function(data_input) {
  # 1. Primary pathway: PMDS -> GDI -> LDL-C
  message("Running causal mediation for GDI (Bootstrap = 1000)...")
  med_gdi <- run_bootstrap_mediation(data_input, mediator_name = "GDI_Score", n_sims = 1000)
  
  # 2. Guild 1 pathway: PMDS -> Guild 1 -> LDL-C
  message("Running causal mediation for Guild 1 (Bootstrap = 1000)...")
  med_g1  <- run_bootstrap_mediation(data_input, mediator_name = "Guild_1", n_sims = 1000)
  
  final_summary <- rbind(med_gdi$summary_table, med_g1$summary_table)
  print(final_summary)
  
  return(final_summary)
}
