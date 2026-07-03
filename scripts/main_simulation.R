## Main simulation script: run all three methods on a single dataset

library(MASS)
library(lme4)
library(glmmTMB)
library(parallel)
library(doParallel)
library(foreach)

source("R/generate_data.R")
source("R/method_wrappers.R")

# Toggle simulation mode: TRUE = fast Test Mode (2 sims), FALSE = Production Mode (100 sims)
run_test_mode <- TRUE

if (run_test_mode) {
  N_SIMS <- 1
  n <- 500
  p <- 8
  n_signals <- 3
  n_groups <- 20
  grV <- 1
  slope <- 1
  sigma_eps <- 1
} else {
  N_SIMS <- 3
  n <- 500
  p <- 8
  n_signals <- 3
  n_groups <- 20
  grV <- 1
  slope <- 1
  sigma_eps <- 1
}

covX <- diag(p)
covU <- diag(slope + 1)

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
  if (is.null(beta_est)) {
    return(out)
  }

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


# Define a single simular function
# Revised iteration function with timing measurements
run_one_iteration_timed <- function(i) {
  
  sim <- generateData(
    beta = beta_true,
    n = n,
    covX = covX,
    n_groups = n_groups,
    grV = grV,
    slope = slope,
    covU = covU,
    sigma_eps = sigma_eps,
    n_candidate_slopes = 3
  )

  data <- sim$data
  data$id <- data$G1

  full_formula <- as.formula(paste(
    "y ~", paste(paste0("x", 1:p), collapse = " + "),
    "+ (1 + x1 | id)"
  ))

  # Helper to record timing and write to file
  record_time <- function(method_name, elapsed) {
    timing_line <- sprintf("Iteration %d, %s, %f seconds", i, method_name, elapsed)
    timing_path <- file.path("results", "timing.txt")
    if (!dir.exists(dirname(timing_path))) dir.create(dirname(timing_path), recursive = TRUE)
    cat(timing_line, "\n", file = timing_path, append = TRUE)
  }

  # glmmTMB (kept on CV)
  t_start <- proc.time()[[3]]
  glmmTMB_fit <- tryCatch({
    run_glmmTMB(data = data, formula = full_formula)
  }, error = function(e) {
    message("glmmTMB failed on iteration ", i, ": ", e$message)
    list(beta = rep(NA_real_, p + 1))
  })
  elapsed <- proc.time()[[3]] - t_start
  record_time("glmmTMB", elapsed)
  beta_glmmTMB <- extract_coefs(glmmTMB_fit$beta, coef_names)

  # Signal-collapse diagnostic: TRUE if any true signal fixed effect (x1, x2,
  # x3) came back below the same 1e-3 hard-threshold fitModel.R applies
  # internally (see R/fitModel.R's `threshold` argument); NA if the fit
  # itself failed (beta_glmmTMB all NA), so a failed fit is never
  # misreported as "no collapse".
  signal_coefs_glmmTMB <- beta_glmmTMB[c("x1", "x2", "x3")]
  glmmTMB_signal_collapse <- if (any(is.na(signal_coefs_glmmTMB))) {
    NA
  } else {
    any(abs(signal_coefs_glmmTMB) < 1e-3)
  }

  # Lin method (BIC selection)
  t_start <- proc.time()[[3]]
  lin_fit <- tryCatch({
    run_lin_method(data = data, random_slope_vars = c("x1", "x2", "x3"))
  }, error = function(e) {
    message("Lin failed on iteration ", i, ": ", e$message)
    list(beta = rep(NA_real_, p + 1))
  })
  elapsed <- proc.time()[[3]] - t_start
  record_time("Lin", elapsed)
  beta_lin <- extract_coefs(lin_fit$beta, coef_names)

  # Bondell method (BIC selection)
  t_start <- proc.time()[[3]]
  bondell_fit <- tryCatch({
    run_bondell_method(data = data, random_slope_vars = c("x1", "x2", "x3"))
  }, error = function(e) {
    message("Bondell failed on iteration ", i, ": ", e$message)
    list(beta = rep(NA_real_, p + 1))
  })
  elapsed <- proc.time()[[3]] - t_start
  record_time("Bondell", elapsed)
  beta_bondell <- extract_coefs(bondell_fit$beta, coef_names)

  # 3. Calculate metrics (unchanged from original)
  metrics_glmmTMB <- calc_metrics(beta_glmmTMB, beta_true)
  metrics_lin <- calc_metrics(beta_lin, beta_true)
  metrics_bondell <- calc_metrics(beta_bondell, beta_true)

  # True random effects matrix (gamma_true, n_groups x 4: Intercept, x1, x2, x3),
  # generated natively by generateData() -- see R/generate_data.R's u_mat/tau_true.
  gamma_true <- sim$u_mat

  # Lin random effects metrics
  if (!is.null(lin_fit$D) && !is.null(lin_fit$bhat) && length(lin_fit$D) > 1) {
    est_slopes_lin <- diag(lin_fit$D)[2:4]
    active_true_RE <- c(TRUE, FALSE, FALSE)
    active_est_RE <- est_slopes_lin > 1e-4

    TP_re <- sum(active_true_RE & active_est_RE)
    FP_re <- sum(!active_true_RE & active_est_RE)
    FN_re <- sum(active_true_RE & !active_est_RE)
    TN_re <- sum(!active_true_RE & !active_est_RE)

    Random_TPR_lin <- ifelse((TP_re + FN_re) > 0, TP_re / (TP_re + FN_re), NA)
    Random_FPR_lin <- ifelse((FP_re + TN_re) > 0, FP_re / (FP_re + TN_re), NA)
    gamma_mspe_lin <- mean((lin_fit$bhat - gamma_true)^2)
    # Raw variance-component estimates for the noise slopes (x2, x3), not
    # just the thresholded active/inactive classification above -- lets us
    # report the magnitude of Lin's occasional non-exact-zero residuals.
    re_var_x2_lin <- unname(est_slopes_lin[2])
    re_var_x3_lin <- unname(est_slopes_lin[3])
  } else {
    Random_TPR_lin <- NA
    Random_FPR_lin <- NA
    gamma_mspe_lin <- NA
    re_var_x2_lin <- NA
    re_var_x3_lin <- NA
  }

  # Bondell random effects metrics
  if (!is.null(bondell_fit$full_result$stddev) && !is.null(bondell_fit$bhat)) {
    est_slopes_bondell <- bondell_fit$full_result$stddev[2:4]
    active_true_RE <- c(TRUE, FALSE, FALSE)
    active_est_RE <- est_slopes_bondell > 1e-4

    TP_re <- sum(active_true_RE & active_est_RE)
    FP_re <- sum(!active_true_RE & active_est_RE)
    FN_re <- sum(active_true_RE & !active_est_RE)
    TN_re <- sum(!active_true_RE & !active_est_RE)

    Random_TPR_bondell <- ifelse((TP_re + FN_re) > 0, TP_re / (TP_re + FN_re), NA)
    Random_FPR_bondell <- ifelse((FP_re + TN_re) > 0, FP_re / (FP_re + TN_re), NA)
    gamma_mspe_bondell <- mean((bondell_fit$bhat - gamma_true)^2)
    # Raw variance-component estimates for x2/x3, squared from stddev to
    # variance so they're on the same scale as Lin's diag(D) above.
    re_var_x2_bondell <- unname(est_slopes_bondell[2])^2
    re_var_x3_bondell <- unname(est_slopes_bondell[3])^2
  } else {
    Random_TPR_bondell <- NA
    Random_FPR_bondell <- NA
    gamma_mspe_bondell <- NA
    re_var_x2_bondell <- NA
    re_var_x3_bondell <- NA
  }

  # glmmTMB random effects metrics (prediction only)
  model_obj <- glmmTMB_fit$full_result$model$model
  if (!is.null(model_obj) && inherits(model_obj, "glmmTMB")) {
    ranef_list <- ranef(model_obj)$cond$id
    ranef_mat <- as.matrix(ranef_list)
    gamma_hat_glmmTMB <- cbind(ranef_mat, 0, 0)
    gamma_mspe_glmmTMB <- mean((gamma_hat_glmmTMB - gamma_true)^2)
  } else {
    gamma_mspe_glmmTMB <- NA
  }
  Random_TPR_glmmTMB <- NA
  Random_FPR_glmmTMB <- NA
  # glmmTMB is fit with the true random-effects formula (1 + x1 | id) --
  # it never estimates x2/x3 random slopes at all, so these are genuinely
  # not applicable (NA), not zero.
  re_var_x2_glmmTMB <- NA
  re_var_x3_glmmTMB <- NA

  # Optimal lambda tracking (unchanged)
  if (!is.null(lin_fit$lambda_re)) {
    opt_lambda_fixed_lin <- lin_fit$lambda_fs
    opt_lambda_random_lin <- lin_fit$lambda_re
  } else {
    opt_lambda_fixed_lin <- NA
    opt_lambda_random_lin <- NA
  }

  if (!is.null(bondell_fit$full_result$t.frac)) {
    opt_lambda_fixed_bondell <- bondell_fit$full_result$t.frac
    opt_lambda_random_bondell <- bondell_fit$full_result$t.frac
  } else {
    opt_lambda_fixed_bondell <- NA
    opt_lambda_random_bondell <- NA
  }

  if (!is.null(glmmTMB_fit$full_result$model$lambda)) {
    opt_lambda_fixed_glmmTMB <- glmmTMB_fit$full_result$model$lambda
  } else {
    opt_lambda_fixed_glmmTMB <- NA
  }
  opt_lambda_random_glmmTMB <- NA

  # Return results for this iteration (same structure, plus the diagnostic
  # columns above)
  result_df <- data.frame(
    Iteration = i,
    Method = c("glmmTMB", "Lin", "Bondell"),
    TPR = c(metrics_glmmTMB["TPR"], metrics_lin["TPR"], metrics_bondell["TPR"]),
    FPR = c(metrics_glmmTMB["FPR"], metrics_lin["FPR"], metrics_bondell["FPR"]),
    MSE = c(metrics_glmmTMB["MSE"], metrics_lin["MSE"], metrics_bondell["MSE"]),
    Random_TPR = c(Random_TPR_glmmTMB, Random_TPR_lin, Random_TPR_bondell),
    Random_FPR = c(Random_FPR_glmmTMB, Random_FPR_lin, Random_FPR_bondell),
    gamma_mspe = c(gamma_mspe_glmmTMB, gamma_mspe_lin, gamma_mspe_bondell),
    opt_lambda_fixed = c(opt_lambda_fixed_glmmTMB, opt_lambda_fixed_lin, opt_lambda_fixed_bondell),
    opt_lambda_random = c(opt_lambda_random_glmmTMB, opt_lambda_random_lin, opt_lambda_random_bondell),
    re_var_x2 = c(re_var_x2_glmmTMB, re_var_x2_lin, re_var_x2_bondell),
    re_var_x3 = c(re_var_x3_glmmTMB, re_var_x3_lin, re_var_x3_bondell),
    glmmTMB_signal_collapse = c(glmmTMB_signal_collapse, NA, NA),
    stringsAsFactors = FALSE
  )

  # Incremental write: each iteration writes its own file the moment it
  # finishes, so a crash partway through a long run doesn't lose already-
  # completed iterations. Each iteration gets its own file (not a shared
  # appended file), which avoids any risk of concurrent-write corruption
  # when multiple %dopar% workers finish at nearly the same time.
  partial_dir <- file.path("results", "partial")
  if (!dir.exists(partial_dir)) dir.create(partial_dir, recursive = TRUE)
  write.csv(result_df, file = file.path(partial_dir, sprintf("iteration_%03d.csv", i)),
            row.names = FALSE)

  result_df
}

# Fresh partial-results directory for this run. The incremental per-iteration
# writes inside run_one_iteration_timed() are the crash-resilience mechanism
# for long runs; clear stale files from any previous run first so they can't
# be mistaken for part of this run's output.
partial_dir <- file.path("results", "partial")
unlink(partial_dir, recursive = TRUE)
dir.create(partial_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# Setup Parallel Processing Cluster
# -----------------------------------------------------------------------------
# Only run in parallel if N_SIMS > 1 to avoid parallel startup overhead and cluster conflicts in single-run tests
use_parallel <- N_SIMS > 1

# Prevent conflicts by closing any existing parallel cluster object if it exists
if (exists("cl") && inherits(cl, "cluster")) {
  try(parallel::stopCluster(cl), silent = TRUE)
}

if (use_parallel) {
  num_cores <- parallel::detectCores() - 1L
  if (num_cores < 1L) num_cores <- 1L
  cat(sprintf("Running in parallel using %d cores...\n", num_cores))

  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  parallel::clusterSetRNGStream(cl, iseed = 12345)
  parallel::clusterEvalQ(cl, {
    source("R/generate_data.R")
    source("R/method_wrappers.R")
  })
  all_results <- foreach(
    i = 1:N_SIMS,
    .combine = rbind,
    .packages = c("MASS", "lme4", "glmmTMB", "quadprog")
  ) %dopar% {
    run_one_iteration_timed(i)
  }
  parallel::stopCluster(cl)
} else {
  cat("Running simulation sequentially (N_SIMS = 1 or parallel disabled)...\n")
  set.seed(12345) # For reproducible sequential runs

  # Run sequentially
  all_results <- foreach(i = 1:N_SIMS, .combine = rbind) %do% {
    run_one_iteration_timed(i)
  }
}
write.csv(all_results, file = "data/simulation_results.csv", row.names = FALSE)
summary_results <- aggregate(cbind(TPR, FPR, MSE, Random_TPR, Random_FPR, gamma_mspe, opt_lambda_fixed, opt_lambda_random, glmmTMB_signal_collapse) ~ Method, data = all_results, FUN = function(x) mean(x, na.rm = TRUE), na.action = na.pass)

cat("\n--- Simulation Summary (Average across", N_SIMS, "iterations) ---\n")
summary_results
