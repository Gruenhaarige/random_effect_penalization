library(MASS)
library(lme4)
library(Matrix)

# beta:       wahrer Parametervektor
# n:          Anzahl Beobachtungen
# covX:       Kovarianzmatrix der Kovariablen
# n_groups:   Anzahl der Gruppen/Cluster für die Random effects  
# grV:        Anzahl der Gruppierungsvariablen
# slope:      Anzahl der random slopes
# covU:       Kovarianzmatrix der random effects
# sigma_eps:  sd des Störterms

generateData <- function(beta, n, covX, n_groups, grV, slope, covU, sigma_eps){
  
  signal <- which(beta[-1] != 0)
  noise <- which(beta[-1] == 0)
  
  p1 <- length(signal)
  p2 <- length(noise)
  p <- p1 + p2
  
  ##############################################################################
  # Gruppierungsvariablen
  
  grV_basic <- rep(seq_len(n_groups), length.out = n)
  
  grV_seq <- seq_len(grV)
  
  names_grV_seq <- paste0("G", seq_len(grV))
  groups_list <- lapply(grV_seq, function(x){
    return(sample(grV_basic))
  })
  names(groups_list) <- names_grV_seq
  groups <- do.call(cbind.data.frame, groups_list)
  
  ##############################################################################
  # Kovariablen
  
  X <- mvrnorm(n, numeric(p), covX)
  X1 <- cbind(rep(1, n), X)
  
  colnames(X) <- paste0("x", seq_len(p))
  
  ##############################################################################
  # Random effects
  
  if(slope > 0){
    slope_var_names <- paste0("x", sort(sample(p, slope)))
    slope_vars <- X[,slope_var_names, drop = FALSE]
  } else{
    slope_var_names <- NULL
    slope_vars <- NULL
  }
  
  makeZpart <- function(group, group_var, slope_vars){
    if(is.null(slope_vars)){
      nsv <- 0
    } else{
      nsv <- NCOL(slope_vars)
    }
    
    Z_intercept <- as.numeric(group_var == group)
    Z_slopes <- lapply(seq_len(nsv), function(isv){
      ifelse(group_var == group, slope_vars[,isv], 0)
    })
    unname(cbind(Z_intercept, do.call(cbind, Z_slopes)))
  }
  
  Z_list <- vector(mode = "list", length = grV)
  u_list <- vector(mode = "list", length = grV)
  
  for(gvar in seq_len(grV)){
    Zp_list <- lapply(seq_len(n_groups), function(g){
      makeZpart(group = g, group_var = groups[[gvar]], slope_vars = slope_vars)
    })
    Z_list[[gvar]] <- do.call(cbind, Zp_list)
    up <- mvrnorm(n_groups, numeric(slope + 1), covU)
    u_list[[gvar]] <- as.vector(t(up))
  }
  
  Z <- do.call(cbind, Z_list)
  u <- do.call(c, u_list)
  
  ##############################################################################
  
  y <- X1 %*% beta + Z %*% u + rnorm(n, sd = sigma_eps)
  list(data = cbind.data.frame(y, X, groups),
       X1 = X1,
       Z = Z,
       u = u,
       slope_var_names = slope_var_names,
       slope_vars = slope_vars)
}

lasso_fixed_effects_objective <- function(beta_full, Y, X_fixed, V_inv, lambda) {
  beta_slopes <- beta_full[-1]
  lasso_penalty <- lambda * sum(abs(beta_slopes))
  residuals <- Y - (X_fixed %*% beta_full)
  GLS_Loss <- t(residuals) %*% V_inv %*% residuals
  objective_value <- GLS_Loss[1,1] + lasso_penalty
  return(objective_value)
}

calculate_V_inv_simple <- function(data, group_factor, sigma_u0, sigma_eps) {
  n <- nrow(data)
  Z_matrix <- model.matrix(~ 0 + as.factor(data[[group_factor]]))
  V <- as.matrix(sigma_u0^2 * (Z_matrix %*% t(Z_matrix)) + sigma_eps^2 * diag(n))
  V_inv <- solve(V)
  return(V_inv)
}

calculate_V_inv_complex <- function(data, model_fit, g2) {
  vc <- VarCorr(model_fit)
  sigma_u_g1 <- vc$G1[1:dim(vc$G1)[1], 1:dim(vc$G1)[2]]
  if (g2 == TRUE){
    sigma_u_g2 <- vc$G2[1:dim(vc$G2)[1], 1:dim(vc$G2)[2]]
  }
  
  sigma_eps_sq <- attr(vc, "sc")^2 
  
  n <- nrow(data)
  
  Z_total <- getME(model_fit, "Z")
  
  n_g1_groups <- length(unique(data$G1))
  Sigma_u_g1_block <- bdiag(rep(list(sigma_u_g1), n_g1_groups))
  
  if (g2 == TRUE){
    n_g2_groups <- length(unique(data$G2))
    Sigma_u_g2_block <- bdiag(rep(list(sigma_u_g2), n_g2_groups))
    Sigma_u_total <- bdiag(Sigma_u_g1_block, Sigma_u_g2_block)
  } else {
    Sigma_u_total <- bdiag(Sigma_u_g1_block)
  }
  
  V <- as.matrix(Z_total %*% Sigma_u_total %*% t(Z_total) + sigma_eps_sq * Diagonal(n))
  V_inv <- solve(V)
  
  return(V_inv)
}
