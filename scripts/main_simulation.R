## Main simulation script: run all three methods on a single dataset

library(MASS)
library(lme4)
library(glmmTMB)

source("R/generate_data.R")
source("R/method_wrappers.R")

# Simulation parameters

{
  N_SIMS <- 2
  
  n <- 500
  p <- 8
  n_signals <- 3    
  n_groups <- 20
  grV <- 1        
  slope <- 1       
  sigma_eps <- 1
  
  covX <- diag(p)
  covU <- diag(slope + 1)
}

# Generate true fixed effects
{
  beta_true <- numeric(p + 1)
  signal_idx <- seq_len(n_signals)
  beta_true[-1][signal_idx] <- 0.5 * (-1)^(signal_idx + 1L)
  names(beta_true) <- c("(Intercept)", paste0("x", 1:p))
  coef_names <- names(beta_true)
}

# Helper function to extract coefficients
extract_coefs <- function(beta_est, coef_names) {
  out <- rep(NA_real_, length(coef_names))
  names(out) <- coef_names
  if (is.null(beta_est)) return(out)
  
  if (is.list(beta_est) && "fixed" %in% names(beta_est)) {
    beta_est <- beta_est$fixed
  } else if (is.list(beta_est) && !is.data.frame(beta_est)) {
    beta_est <- unlist(beta_est)
  }
  
  if (is.null(names(beta_est))) {
    len <- min(length(beta_est), length(coef_names))
    out[seq_len(len)] <- beta_est[seq_len(len)]
  } else {
    common <- intersect(names(beta_est), coef_names)
    out[common] <- beta_est[common]
  }
  out
}

# Helper function to calculate metrics
calc_metrics <- function(est, truth) {
  # Exclude intercept for TPR/FPR
  est_vars <- est[-1]
  truth_vars <- truth[-1]
  
  mse <- mean((est - truth)^2, na.rm = TRUE)
  
  active_true <- truth_vars != 0
  active_est <- !is.na(est_vars) & abs(est_vars) > 1e-4 
  
  TP <- sum(active_true & active_est)
  FP <- sum(!active_true & active_est)
  FN <- sum(active_true & !active_est)
  TN <- sum(!active_true & !active_est)
  
  TPR <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
  FPR <- ifelse((FP + TN) > 0, FP / (FP + TN), NA)
  
  return(c(TPR = TPR, FPR = FPR, MSE = mse))
}


results_list <- list()

# Simulation Loop
for (i in 1:N_SIMS) {
  cat(sprintf("Iteration %d / %d\n", i, N_SIMS))
  
  # 1. Generate data
  sim <- generateData(beta = beta_true,
                      n = n,
                      covX = covX,
                      n_groups = n_groups,
                      grV = grV,
                      slope = slope,
                      covU = covU,
                      sigma_eps = sigma_eps)
  
  data <- sim$data
  data$id <- data$G1
  
  full_formula <- as.formula(paste("y ~", paste(paste0("x", 1:p), collapse = " + "),
                                   "+ (1 + x1 | id)"))
  
  # 2. Apply methods
  
  # glmmTMB
  glmmTMB_fit <- tryCatch({ run_glmmTMB(data = data, formula = full_formula) }, error = function(e) list(beta = rep(NA, p+1)))
  beta_glmmTMB <- extract_coefs(glmmTMB_fit$beta, coef_names)
  
  # Lin method
  lin_fit <- tryCatch({ run_lin_method(data = data) }, error = function(e) list(beta = rep(NA, p+1)))
  beta_lin <- extract_coefs(lin_fit$beta, coef_names)
  
  # Bondell method
  bondell_fit <- tryCatch({ run_bondell_method(data = data, random_slope_var = "x1") }, error = function(e) list(beta = rep(NA, p+1)))
  beta_bondell <- extract_coefs(bondell_fit$beta, coef_names)
  
  # 3. Calculate metrics
  metrics_glmmTMB <- calc_metrics(beta_glmmTMB, beta_true)
  metrics_lin <- calc_metrics(beta_lin, beta_true)
  metrics_bondell <- calc_metrics(beta_bondell, beta_true)
  
  # 4. Save results
  results_list[[i]] <- data.frame(
    Iteration = i,
    Method = c("glmmTMB", "Lin", "Bondell"),
    TPR = c(metrics_glmmTMB["TPR"], metrics_lin["TPR"], metrics_bondell["TPR"]),
    FPR = c(metrics_glmmTMB["FPR"], metrics_lin["FPR"], metrics_bondell["FPR"]),
    MSE = c(metrics_glmmTMB["MSE"], metrics_lin["MSE"], metrics_bondell["MSE"])
  )
}

all_results <- do.call(rbind, results_list)
write.csv(all_results, file = "data/simulation_results.csv", row.names = FALSE)
summary_results <- aggregate(cbind(TPR, FPR, MSE) ~ Method, data = all_results, FUN = function(x) mean(x, na.rm = TRUE))

cat("\n--- Simulation Summary (Average across", N_SIMS, "iterations) ---\n")
summary_results
