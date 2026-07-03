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
  measure <- "AIC"
  n_lambda <- 100
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
                           lambda_re = 5,
                           lambda_fs = 10,
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
    matrix(1, nrow = nrow(df), ncol = 1)
  })
  
  y_combined <- do.call(rbind, y_list)
  x_combined <- do.call(rbind, x_list)
  z_combined <- do.call(rbind, z_list)
  
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
  
  ## STAGE 1: Pen.fs (random effects selection)
  pe1 <- tryCatch({
    Pen.fs(lambda = lambda_re, 
           x = x_list, 
           y = y_list, 
           z = z_list,
           D.init = D.init, 
           sig.init = sig.init, 
           eps = eps_reml)
  }, error = function(e) {
    stop("Pen.fs (Stage 1) failed: ", e$message)
  })
  
  ## STAGE 2: pco (fixed effects selection)
  final_fit <- tryCatch({
    pco(x = x_list, 
        y = y_list, 
        z = z_list,
        beta.init = pe1$beta, 
        D.init = pe1$D, 
        sig.init = pe1$sig,
        lambda = lambda_fs, 
        eps = eps_pco)
  }, error = function(e) {
    stop("pco (Stage 2) failed: ", e$message)
  })
  
  beta_est <- final_fit$beta
  
  if (length(beta_est) == length(x_cols) + 1) {
    names(beta_est) <- c("(Intercept)", x_cols)
  }
  
  return(list(
    beta = beta_est,
    D = final_fit$D,
    sig = final_fit$sig,
    stage1_result = pe1
  ))
}

## -----------------------------------------------------------------------------
## Bondell's Pen.LME method
## -----------------------------------------------------------------------------
## Uses Pen.LME

run_bondell_method <- function(data,
                                random_slope_var = "x1",
                                t.fracs = seq(1, 0.2, -0.2),
                                ...) {
  
  x_cols <- grep("^x[0-9]+$", colnames(data), value = TRUE)
  
  subject <- as.factor(data$id)
  
  fixed_formula <- as.formula(paste("y ~", paste(x_cols, collapse = " + ")))
  
  y_vec <- data$y
  X_mat <- model.matrix(fixed_formula, data = data)
  
  ## Construct Z matrix for random effects
  ## Pen.LME automatically adds an intercept column, so we only need the slope variables
  if (!is.null(random_slope_var) && random_slope_var %in% colnames(data)) {
    Z_mat <- as.matrix(data[[random_slope_var]], ncol = 1)
  } else {
    ## Random intercept only: pass empty matrix (intercept will be added by Pen.LME)
    Z_mat <- matrix(nrow = nrow(data), ncol = 0)
  }
  
  result <- Pen.LME(y = y_vec,
                    X = X_mat,
                    Z = Z_mat,
                    subject = subject,
                    t.fracs = t.fracs)
  
  if (is.list(result) && "beta" %in% names(result)) {
    beta_est <- result$beta
  } else if (is.vector(result)) {
    beta_est <- result
  } else {
    beta_est <- result$coefficients
  }
  
  if (length(beta_est) == ncol(X_mat)) {
    names(beta_est) <- colnames(X_mat)
  }
  
  return(list(
    beta = beta_est,
    full_result = result
  ))
}

