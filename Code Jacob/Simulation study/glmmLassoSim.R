glmmLassoSim <- function(n,
                         p1,
                         p2,
                         intercept = 0,
                         beta_part = c(0.5, -0.5),
                         covX_strength = 0,
                         sigma_eps = 1,
                         slope = 0L,
                         n_grouping_var = 1L,
                         n_groups = 5L,
                         covU_strength = 0.3,
                         method = "glmmTMB",
                         common_scale = TRUE,
                         lambda_max_start = 10,
                         n_lambda = 100L,
                         const = 1e-5,
                         location = 0,
                         threshold = 1e-3,
                         max_iter = 100000,
                         max_run_lambda_search = 10L,
                         message_value = 10^floor(log(0.2 * n_lambda, 10))){
  
  p <- p1 + p2
  
  ##############################################################################
  # internal functions:

  is.natural <- function(x){
    is.numeric(x) && all(x >= 0) && all(x %% 1 == 0)
  }

  getCovMat <- function(d, strength){
    out <- matrix(strength, d, d)
    diag(out) <- 1
    out
  }

  ##############################################################################
  # input checks:

  stopifnot(is.natural(n),
            is.natural(p1),
            is.natural(p2),
            is.numeric(intercept),
            is.numeric(beta_part),
            is.numeric(covX_strength),
            covX_strength >= 0,
            is.numeric(sigma_eps),
            sigma_eps > 0,
            length(slope) == 1L,
            is.natural(slope),
            slope <= p,
            is.natural(n_grouping_var),
            is.natural(n_groups),
            is.numeric(covU_strength),
            covU_strength >= 0,
            method %in% c("blme", "glmmTMB", "lme4"),
            is.logical(common_scale),
            length(lambda_max_start) == 1L,
            is.numeric(lambda_max_start),
            lambda_max_start > 0,
            is.natural(n_lambda)
  )

  ##############################################################################
  # data generation:

  if(length(beta_part) > p1){
    warning("Cannot use the full beta_part!")
  }

  true_beta <- c(intercept,
                 rep(beta_part, length.out = p1),
                 numeric(p2))

  q <- slope + 1L

  dataList <- generateData(beta = true_beta,
                           n = n,
                           covX = getCovMat(p, covX_strength),
                           n_groups = n_groups,
                           grV = n_grouping_var,
                           slope = slope,
                           covU = getCovMat(q, covU_strength),
                           sigma_eps = sigma_eps)
  
  simData <- dataList[["data"]]
  X1 <- dataList[["X1"]]
  Z <- dataList[["Z"]]
  u <- dataList[["u"]]

  ##############################################################################
  # building formula:

  comb_slope <- expand.grid(Xvar = dataList$slope_var_names,
                            Gvar = paste0("G", seq_len(n_grouping_var))
  )
  ncs <- nrow(comb_slope)

  covariates <- paste0("x", seq_len(p))

  prepare_formula <- c(
    covariates,
    paste0("(1|G", seq_len(n_grouping_var), ")")
  )

  if(slope > 0){
    prepare_formula <- c(
      prepare_formula,
      sapply(seq_len(ncs), function(k){
        paste0("(", comb_slope$Xvar[k], "|", comb_slope$Gvar[k], ")")
      })
    )
  }

  formula <- formula(paste("y ~", paste(prepare_formula, collapse = " + ")))
  
  ##############################################################################
  
  if(method %in% c("blme", "glmmTMB")){
    model <- fitGlmmLasso(formula, data = simData, method = method,
                          lambda_max_start = lambda_max_start,
                          n_lambda = n_lambda, const = const,
                          location = location, threshold = threshold,
                          max_iter = max_iter,
                          lambda_max_runs = max_run_lambda_search,
                          message_value = message_value)
    
    return(list(
      overview = model$overview,
      true_beta = true_beta,
      param = model$param,
      slope_var_names = dataList$slope_var_names
    ))
  } else if(method == "lme4"){
    model <- fitModel(formula, data = simData, method = "lme4",
                      max_iter = max_iter)
    
    return(list(
      overview = NULL,
      true_beta = true_beta,
      param = model$beta,
      slope_var_names = dataList$slope_var_names
    ))
  } else{
    stop("Unknown method!")
  }
}