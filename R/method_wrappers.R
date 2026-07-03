# Method wrapper functions for simulation comparisons
#
# run_glmmTMB(): glmmTMB-based method with lambda search
# run_lin_method(): Lin's two-stage method (Pen.fs + pco)
# run_bondell_method(): Bondell's Pen.LME method

library(MASS)
library(lme4)
library(glmmTMB)
library(Matrix)

source("R/fitModel.R")
source("R/fitGlmmLasso.R")
source("R/method_lin.R")
source("R/method_bondell.R")

# -----------------------------------------------------------------------------
# glmmTMB-based method
# -----------------------------------------------------------------------------

run_glmmTMB <- function(data, formula) {
  # response_name <- all.vars(formula)[1]
  # varnames <- attr(terms(lme4::nobars(formula)), "term.labels")
  
  method <- "glmmTMB"
  measure <- "CV"
  n_lambda <- 10
  lambda_max_adjust <- TRUE
  lambda_max_start <- 100
  lambda_max_runs <- 10
  
  result <- fitGlmmLasso(
    formula = formula,
    data = data,
    method = method,
    measure = measure,
    n_lambda = n_lambda,
    lambda_max_adjust = lambda_max_adjust,
    lambda_max_start = lambda_max_start,
    lambda_max_runs = lambda_max_runs
  )
  
  if (!is.na(result$param[1])) {
    beta_est <- result$param
  } else {
    beta_est <- NA
    warning("glmmTMB method failed to converge.")
  }
  
  return(list(
    beta = beta_est,
    full_result = result
  ))
}

## -----------------------------------------------------------------------------
## Lin's two-stage method
## -----------------------------------------------------------------------------
## Stage 1: Pen.fs (penalized random effects selection)
## Stage 2: pco (penalized fixed effects selection)

run_lin_method <- function(data, 
                           random_slope_vars = c("x1", "x2", "x3"),
                           eps_reml = 1e-5,
                           eps_pco = 1e-5) {
  
  data_list <- split(data, data$id)
  n_subjects <- length(data_list)
  
  x_cols <- grep("^x[0-9]+$", colnames(data), value = TRUE)
  
  y_list <- lapply(data_list, function(df) as.matrix(df$y))
  x_list <- lapply(data_list, function(df) {
    X_mat <- as.matrix(df[, x_cols, drop = FALSE])
    X_with_int <- cbind(1, X_mat)
    colnames(X_with_int) <- c("(Intercept)", x_cols)
    X_with_int
  })
  
  z_list <- lapply(data_list, function(df) {
    Z_mat <- cbind(1, as.matrix(df[, random_slope_vars, drop = FALSE]))
    colnames(Z_mat) <- c("(Intercept)", random_slope_vars)
    Z_mat
  })
  
  y_combined <- do.call(rbind, y_list)
  x_combined <- do.call(rbind, x_list)
  z_combined <- do.call(rbind, z_list)
  n_total <- nrow(y_combined)

  subject_ids <- rep(names(data_list), times = sapply(data_list, nrow))
  subject <- as.factor(subject_ids)

  lmer_data <- data.frame(
    y_response = y_combined,
    x_matrix = I(x_combined),
    z_matrix = I(z_combined),
    subject_id = subject
  )

  init_fit <- lmer(y_response ~ x_matrix - 1 + (0 + z_matrix | subject_id),
                   data = lmer_data, REML = TRUE)

  hh <- VarCorr(init_fit)
  sig.init <- (attributes(hh)$sc)^2
  D.init <- as.matrix(hh[[1]])

  # Grids traversed largest -> smallest, with warm starts at every step
  # (De.start/sig.start/beta.start below carry the previous grid point's
  # converged fit into the next one).
  lambda_re_grid <- c(40, 20, 10, 5, 2.5, 1, 0.5, 0.1)
  lambda_fs_grid <- c(80, 40, 20, 10, 5, 2.5, 1, 0.5, 0.1)

  # -------------------------------------------------------------------------
  # Stage 1: BIC-select lambda_re over the REML variance-component fit.
  # D is diagonal by assumption (Lin's method estimates no random-effect
  # covariances; contrast with Bondell's full Cholesky-based D), so df counts
  # nonzero variance components on the diagonal, not lower-triangle entries.
  # BIC = -2*REML_logLik + log(m)*df, m = number of clusters.
  # Each grid point's fit is kept in stage1_fits so the final model is
  # exactly the fit that produced the selected BIC value (no separate,
  # non-warm-started refit).
  # -------------------------------------------------------------------------
  De_curr <- D.init
  sig_curr <- sig.init
  stage1_fits <- vector("list", length(lambda_re_grid))
  bic_re <- numeric(length(lambda_re_grid))

  for (k in seq_along(lambda_re_grid)) {
    fit_k <- Pen.fs(lambda = lambda_re_grid[k], x = x_list, y = y_list, z = z_list,
                     D.init = D.init, sig.init = sig.init, eps = eps_reml,
                     De.start = De_curr, sig.start = sig_curr)
    stage1_fits[[k]] <- fit_k
    De_curr <- fit_k$D
    sig_curr <- fit_k$sig

    df1 <- sum(diag(fit_k$D) != 0) + 1
    ll1 <- reml_loglik_lin(x_list, y_list, z_list, fit_k$beta, fit_k$D, fit_k$sig)
    bic_re[k] <- -2 * ll1 + log(n_subjects) * df1
  }
  if (all(!is.finite(bic_re))) {
    stop("run_lin_method: no finite BIC values across lambda_re grid (Stage 1).")
  }
  bic_re[!is.finite(bic_re)] <- Inf

  k_star <- which.min(bic_re)
  lambda_re_opt <- lambda_re_grid[k_star]
  pe1_opt <- stage1_fits[[k_star]]

  # -------------------------------------------------------------------------
  # Stage 2: BIC-select lambda_fs over the adaptive-lasso fixed-effect fit,
  # holding D/sig fixed at the Stage-1 optimum (pe1_opt).
  # BIC = -2*logLik + log(N)*df, N = total observations.
  # -------------------------------------------------------------------------
  beta_curr <- pe1_opt$beta
  stage2_fits <- vector("list", length(lambda_fs_grid))
  bic_fs <- numeric(length(lambda_fs_grid))

  for (j in seq_along(lambda_fs_grid)) {
    fit_j <- pco(x = x_list, y = y_list, z = z_list,
                 beta.init = pe1_opt$beta, D.init = pe1_opt$D, sig.init = pe1_opt$sig,
                 lambda = lambda_fs_grid[j], eps = eps_pco, beta.start = beta_curr)
    stage2_fits[[j]] <- fit_j
    beta_curr <- fit_j$beta

    df2 <- sum(fit_j$beta != 0)
    ll2 <- marginal_loglik_lin(x_list, y_list, z_list, fit_j$beta, pe1_opt$D, pe1_opt$sig)
    bic_fs[j] <- -2 * ll2 + log(n_total) * df2
  }
  if (all(!is.finite(bic_fs))) {
    stop("run_lin_method: no finite BIC values across lambda_fs grid (Stage 2).")
  }
  bic_fs[!is.finite(bic_fs)] <- Inf

  j_star <- which.min(bic_fs)
  lambda_fs_opt <- lambda_fs_grid[j_star]
  final_fit <- stage2_fits[[j_star]]

  beta_est <- as.vector(final_fit$beta)
  if (length(beta_est) == length(x_cols) + 1) {
    names(beta_est) <- c("(Intercept)", x_cols)
  }
  
  bhat_lin <- matrix(0, nrow = n_subjects, ncol = 4)
  for (i in 1:n_subjects) {
    y_i <- y_list[[i]]
    X_i <- x_list[[i]]
    Z_i <- z_list[[i]]
    V_i <- Z_i %*% final_fit$D %*% t(Z_i) + final_fit$sig * diag(nrow(Z_i))
    V_i_inv <- ginv(V_i)
    bhat_lin[i, ] <- final_fit$D %*% t(Z_i) %*% V_i_inv %*% (y_i - X_i %*% beta_est)
  }
  
  return(list(
    beta = beta_est,
    D = final_fit$D,
    sig = final_fit$sig,
    stage1_result = pe1_opt,
    bhat = bhat_lin,
    lambda_re = lambda_re_opt,
    lambda_fs = lambda_fs_opt,
    bic_re = bic_re,
    bic_fs = bic_fs
  ))
}

## -----------------------------------------------------------------------------
## Bondell's Pen.LME method
## -----------------------------------------------------------------------------
## Uses Pen.LME

run_bondell_method <- function(data,
                                random_slope_vars = c("x1", "x2", "x3"),
                                t.fracs = seq(1, 0.2, -0.2),
                                ...) {
  
  # Crucial step: Pen.LME expects data to be sorted/grouped by subject ID
  data <- data[order(data$id), ]
  
  x_cols <- grep("^x[0-9]+$", colnames(data), value = TRUE)
  
  subject <- as.factor(data$id)
  
  fixed_formula <- as.formula(paste("y ~", paste(x_cols, collapse = " + ")))
  
  y_vec <- data$y
  X_mat <- model.matrix(fixed_formula, data = data)
  
  ## Construct Z matrix for random effects
  ## Pen.LME automatically adds an intercept column, so we only need the slope variables
  if (!is.null(random_slope_vars) && all(random_slope_vars %in% colnames(data))) {
    Z_mat <- as.matrix(data[, random_slope_vars, drop = FALSE])
  } else {
    ## Random intercept only: pass empty matrix (intercept will be added by Pen.LME)
    Z_mat <- matrix(nrow = nrow(data), ncol = 0)
  }
  
  result <- Pen.LME(y = y_vec,
                    X = X_mat,
                    Z = Z_mat,
                    subject = subject,
                    t.fracs = t.fracs)
  
  if (is.list(result) && "fixed" %in% names(result)) {
    beta_est <- result$fixed
  } else if (is.list(result) && "beta" %in% names(result)) {
    beta_est <- result$beta
  } else if (is.vector(result)) {
    beta_est <- result
  } else {
    beta_est <- result$coefficients
  }
  
  beta_est <- as.vector(beta_est)
  
  if (length(beta_est) == ncol(X_mat)) {
    names(beta_est) <- colnames(X_mat)
  }
  
  return(list(
    beta = beta_est,
    full_result = result,
    bhat = result$bhat
  ))
}


# End of method_wrappers.R

NULL

