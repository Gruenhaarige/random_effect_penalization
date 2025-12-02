fitGlmmLasso <- function(formula, data, method, family = gaussian,
                         measurement = "AIC", n_lambda = 100,
                         lambda_max_adjust = TRUE, lambda_max_start = 10,
                         lambda_max_runs = 10, lsl = 5, max_iter = 100000,
                         const = 1e-5, location = 0, threshold = 1e-3, svf = 1,
                         max_attempts = c(10, 12),
                         message_value = 10^floor(log(0.2 * n_lambda, 10))){
  
  # max_attempts: (without, with) varying lambda
  # lsl = lambda_seq_level: bendinglevel of lambda_seq
  # svf = starting value factor: new = svf * est_beta + (1-svf) * intercept_model
  
  ##############################################################################
  # find lambda_max & lambda_seq
  
  response_name <- all.vars(formula)[1]
  varnames <- attr(terms(nobars(formula)), "term.labels")
  default_start <- c(mean(data[[response_name]]), numeric(length(varnames)))
  
  first_fit <- fitModel(formula = formula,
                        data = data,
                        method = method,
                        lambda = lambda_max_start,
                        const = const,
                        location = location,
                        threshold = threshold,
                        max_iter = max_iter,
                        start = default_start,
                        vary_lambda = lambda_max_adjust,
                        max_attempts = ifelse(lambda_max_adjust,
                                              max_attempts[2],
                                              max_attempts[1]))
  
  test_beta <- first_fit$beta
  k <- 0L
  
  if(any(test_beta[-1] != 0)){
    lambda_old <- lambda_max_start
    lambda_new <- first_fit$lambda * 2
    message("lambda = ", lambda_old, " too small. Try lambda = ", lambda_new)
    repeat{
      test_fit <- fitModel(formula = formula,
                           data = data,
                           method = method,
                           lambda = lambda_new,
                           const = const,
                           location = location,
                           threshold = threshold,
                           common_scale = TRUE,
                           max_iter = max_iter,
                           start = test_beta,
                           vary_lambda = lambda_max_adjust,
                           max_attempts = ifelse(lambda_max_adjust,
                                                 max_attempts[2],
                                                 max_attempts[1]))
      k <- k + 1L
      
      if(k == lambda_max_runs){
        lambda_max <- lambda_new
        warning("The specified value for 'lambda_max' may be too small!")
        break
      } else{
        if(any(test_fit$beta[-1] != 0)){
          lambda_old <- lambda_new
          lambda_new <- test_fit$lambda * 2
          test_beta <- test_fit$beta
          message("lambda = ", lambda_old, " too small. Try lambda = ",
                  lambda_new)
          last_success <- lambda_old
        } else{
          lambda_mid <- (lambda_old + lambda_new) / 2
          test_beta <- test_fit$beta
          final_fit <- fitModel(formula = formula,
                                data = data,
                                method = method,
                                lambda = lambda_mid,
                                const = const,
                                location = location,
                                threshold = threshold,
                                common_scale = TRUE,
                                max_iter = max_iter,
                                start = test_beta,
                                vary_lambda = lambda_max_adjust,
                                max_attempts = ifelse(lambda_max_adjust,
                                                      max_attempts[2],
                                                      max_attempts[1]))
          
          if(any(final_fit$beta[-1] != 0)){
            lambda_max <- lambda_new
          } else{
            lambda_max <- lambda_mid
          }
          break
        }
      }
    }
  } else{
    lambda_old <- lambda_max_start
    lambda_new <- first_fit$lambda / 2
    message("lambda = ", lambda_old, " may be too large. Try lambda = ",
            lambda_new)
    repeat{
      test_fit <- fitModel(formula = formula,
                           data = data,
                           method = method,
                           lambda = lambda_new,
                           const = const,
                           location = location,
                           threshold = threshold,
                           common_scale = TRUE,
                           max_iter = max_iter,
                           start = test_beta,
                           vary_lambda = lambda_max_adjust,
                           max_attempts = ifelse(lambda_max_adjust,
                                                 max_attempts[2],
                                                 max_attempts[1]))
      k <- k + 1L
      
      if(k == lambda_max_runs){
        if(any(test_fit$beta[-1] != 0)){
          lambda_max <- lambda_old
          warning("The specified value for 'lambda_max' may be too large!")
        } else{
          lambda_max <- lambda_new
        }
        break
      } else{
        if(any(test_fit$beta[-1] != 0)){
          lambda_mid <- (lambda_old + lambda_new) / 2
          test_beta <- test_fit$beta
          final_fit <- fitModel(formula = formula,
                                data = data,
                                method = method,
                                lambda = lambda_mid,
                                const = const,
                                location = location,
                                threshold = threshold,
                                common_scale = TRUE,
                                max_iter = max_iter,
                                start = test_beta,
                                vary_lambda = lambda_max_adjust,
                                max_attempts = ifelse(lambda_max_adjust,
                                                      max_attempts[2],
                                                      max_attempts[1]))
          
          if(any(final_fit$beta[-1] != 0)){
            lambda_max <- lambda_old
          } else{
            lambda_max <- lambda_mid
          }
          break
        } else{
          lambda_old <- lambda_new
          lambda_new <- test_fit$lambda / 2
          test_beta <- test_fit$beta
          message("lambda = ", lambda_old, " may be too large. Try lambda = ",
                  lambda_new)
          last_success <- lambda_old
        }
      }
    }
  }
  message("Found lambda_max: ", lambda_max)
  
  # lambda_seq <- exp(seq(log(lambda_max), log(1e-3), length.out = n_lambda))
  lambda_seq <- seq(lambda_max^(1 / lsl), 1e-3^(1 / lsl), length.out = n_lambda)
  lambda_seq <- lambda_seq^lsl
  
  ##############################################################################
  # calculate beta paths
  
  nonzero <- logical(n_lambda)
  
  param <- vector(mode = "list", length = n_lambda)
  names(param) <- paste0("l", lambda_seq)
  
  n_nonzero_param <- numeric(n_lambda)
  LL <- numeric(n_lambda)
  AIC <- numeric(n_lambda)
  BIC <- numeric(n_lambda)
  
  fail <- character(n_lambda)
  new_start <- default_start
  
  for(i in seq_len(n_lambda)){
    pen <- lambda_seq[i]
    nfail <- 0L
    
    model <- fitModel(formula = formula,
                      data = data,
                      method = method,
                      lambda = pen,
                      const = const,
                      location = location,
                      threshold = threshold,
                      common_scale = TRUE,
                      max_iter = max_iter,
                      start = new_start,
                      max_attempts = max_attempts[1])
    
    param[[i]] <- model$beta
    new_start <- svf * model$beta + (1 - svf) * default_start
    LL[i] <- model$LL
    n_nonzero_param[i] <- model$nonzero
    nonzero[i] <- any(model$beta[-1] != 0)
    
    if(i %% message_value == 0){
      message(i, " / ", n_lambda)
    }
  }
  
  AIC <- -2 * LL + 2 * n_nonzero_param
  BIC <- -2 * LL + log(nrow(data)) * n_nonzero_param
  
  ##############################################################################
  # optimal model & output
  
  if(measurement == "AIC"){
    lambda_opt_i <- which.min(AIC)
  } else if(measurement == "BIC"){
    lambda_opt_i <- which.min(BIC)
  } else{
    stop("unknown measurement!")
  }
  
  if(length(lambda_opt_i) == 1){
    lambda_opt <- lambda_seq[lambda_opt_i]
    if(lambda_opt_i == 1L){
      final_start <- default_start
    } else{
      final_start <- param[[lambda_opt_i - 1L]]
    }
    
    opt_model <- fitModel(formula = formula,
                          data = data,
                          method = method,
                          lambda = lambda_opt,
                          const = const,
                          location = location,
                          threshold = threshold,
                          common_scale = TRUE,
                          max_iter = max_iter,
                          start = final_start,
                          max_attempts = max_attempts[1])
    
    return(list(
      param = opt_model$beta,
      overview = data.frame(lambda_seq = lambda_seq,
                            logLik = LL,
                            n_nonzero_param = n_nonzero_param,
                            AIC = AIC,
                            BIC = BIC,
                            fail
      ),
      optimal_model = opt_model,
      lambda_max_search_runs = k,
      param_list = param
    ))
    
  } else{
    message("Länge von lambda_opt_i ist ", length(lambda_opt_i), " Objekt:")
    message(lambda_opt_i)
    
    return(list(
      param = NA,
      overview = data.frame(lambda_seq = lambda_seq,
                            logLik = LL,
                            n_nonzero_param = n_nonzero_param,
                            AIC = AIC,
                            BIC = BIC,
                            fail
      ),
      optimal_model = NA,
      lambda_max_search_runs = k,
      param_list = param
    ))
  }
}