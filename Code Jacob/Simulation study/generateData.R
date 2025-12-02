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
  X <- scale(X)
  X1 <- cbind(rep(1, n), X)
  
  colnames(X) <- paste0("x", seq_len(p))
  
  ##############################################################################
  # Random effects
  
  if(slope > 0){
    slope_var_names <- paste0("x", sort(sample(p, slope)))
    slope_vars <- X[,slope_var_names]
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