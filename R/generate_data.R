generateData <- function(beta, n, covX, n_groups, grV, slope, covU, sigma_eps,
                          n_candidate_slopes = slope){

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
  #
  # n_candidate_slopes is the full set of random-slope CANDIDATES considered
  # by the variable-selection methods; it defaults to `slope` so any caller
  # that does not pass it explicitly gets byte-identical behaviour to the
  # old function (zero padding added below). Passing n_candidate_slopes >
  # slope (as main_simulation.R does, with slope = 1, n_candidate_slopes = 3)
  # makes the DGP natively carry the noise candidate slopes (x2, x3) as real
  # regressors in Z with an EXACT ZERO true random-effect coefficient,
  # instead of that structure being hand-reconstructed later. The true
  # slopes (first `slope` of the candidates, matching covU) are drawn with
  # the exact same mvrnorm() call as before; the extra candidate columns are
  # appended afterwards with no additional random draws, so y is generated
  # identically to the old 2-dimensional sampling for the same seed.
  stopifnot(n_candidate_slopes >= slope)

  candidate_slope_names <- if (n_candidate_slopes > 0) paste0("x", seq_len(n_candidate_slopes)) else NULL
  candidate_slope_vars <- if (n_candidate_slopes > 0) X[, candidate_slope_names, drop = FALSE] else NULL

  slope_var_names <- if (slope > 0) paste0("x", seq_len(slope)) else NULL
  slope_vars <- if (slope > 0) X[, slope_var_names, drop = FALSE] else NULL

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

  q_re <- 1 + n_candidate_slopes
  re_names <- c("Intercept", candidate_slope_names)
  n_pad <- n_candidate_slopes - slope

  # True random-effect covariance D (q_re x q_re): covU occupies the
  # intercept + true-slope block (tau0^2 = tau1^2 = 1 by default, unchanged
  # signal strength); the padded candidate slopes have exactly zero true
  # variance and zero covariance with everything else.
  D_true <- matrix(0, nrow = q_re, ncol = q_re, dimnames = list(re_names, re_names))
  D_true[1:(slope + 1), 1:(slope + 1)] <- covU
  tau_true <- diag(D_true)

  Z_list <- vector(mode = "list", length = grV)
  u_list <- vector(mode = "list", length = grV)

  for(gvar in seq_len(grV)){
    Zp_list <- lapply(seq_len(n_groups), function(g){
      makeZpart(group = g, group_var = groups[[gvar]], slope_vars = candidate_slope_vars)
    })
    Z_list[[gvar]] <- do.call(cbind, Zp_list)

    # Sample only the true (nonzero-variance) intercept+slope block, exactly
    # as before, then pad with exact-zero columns for the noise candidates
    # -- this avoids ever drawing from a singular MVN.
    up_true <- mvrnorm(n_groups, numeric(slope + 1), covU)
    up <- cbind(up_true, matrix(0, nrow = n_groups, ncol = n_pad))
    colnames(up) <- re_names
    u_list[[gvar]] <- as.vector(t(up))
  }

  Z <- do.call(cbind, Z_list)
  u <- do.call(c, u_list)

  ##############################################################################

  y <- X1 %*% beta + Z %*% u + rnorm(n, sd = sigma_eps)
  # Create spurious grouping factor by permuting true IDs
  id_fake <- sample(groups$G1)

  # u_mat: n_groups x q_re matrix of the true per-cluster random effects
  # actually used above (columns named by re_names, e.g. Intercept, x1, x2,
  # x3). Only meaningful for the first grouping variable (grV == 1 is the
  # only case used in this simulation study).
  u_mat <- matrix(u_list[[1]], nrow = n_groups, ncol = q_re, byrow = TRUE,
                   dimnames = list(NULL, re_names))

  list(data = cbind.data.frame(y, X, groups, id_fake = id_fake),
       X1 = X1,
       Z = Z,
       u = u,
       u_mat = u_mat,
       tau_true = tau_true,
       D_true = D_true,
       slope_var_names = slope_var_names,
       slope_vars = slope_vars)
}
