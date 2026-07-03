fitGlmmLasso <- function(formula,
                         data,
                         method,
                         family = gaussian,
                         measure = c("CV", "BIC", "AIC"),
                         n_folds = NULL,
                         LL_with_RE = TRUE,
                         n_lambda = 100,
                         lambda_max_adjust = TRUE,
                         lambda_max_start = 100,
                         lambda_max_runs = 15,
                         lsl = 5,
                         max_iter = 100000,
                         const = 1e-5,
                         location = 0,
                         threshold = 1e-3,
                         # lme4_REML = TRUE,
                         svf = 1,
                         max_attempts = c(10, 12),
                         message_value = 10^floor(log(0.2 * n_lambda, 10))){
  
  # max_attempts: (without, with) varying lambda
  # lsl = lambda_seq_level: bendinglevel of lambda_seq
  # svf = starting value factor: new = svf * est_beta + (1-svf) * intercept_model
  
  if(measure[1] == "CV"){
    gvar <- sapply(lme4::findbars(formula), function(x) as.character(x[[3]]))
    ngv <- length(gvar)
    if(ngv > 1L){
      comb_id <- apply(data[, gvar], 1, function(x) paste(x, collapse = "_"))
      unique_id <- unique(comb_id)
      cv_group <- match(comb_id, unique_id)
    } else{
      unique_id <- unique(data[[gvar]])
      cv_group <- match(data[[gvar]], unique_id)
    }
    if(is.null(n_folds)){
      n_folds <- length(unique_id)
    } else{
      stopifnot(length(n_folds) == 1,
                is.numeric(n_folds),
                n_folds > 1L,
                n_folds %% 1 == 0)
      if(n_folds > length(unique_id)){
        stop("Anzahl folds darf nicht größer sein, als die Anzahl der Gruppen!")
      }
    }
  } else if(measure[1] %in% c("BIC", "AIC")){
    n_folds <- 1L
  } else{
    stop("Ungültige Wahl für 'measure'!")
  }
  
  ##############################################################################
  # find lambda_max & lambda_seq
  
  response_name <- all.vars(formula)[1]
  varnames <- attr(terms(nobars(formula)), "term.labels")
  
  
  if(method %in% c("blme", "glmmTMB")){
    default_start <- c(mean(data[[response_name]], na.rm = TRUE),
                       numeric(length(varnames)))
    
    first_fit <- fitModel(formula = formula,
                          data = data,
                          method = method,
                          lambda = lambda_max_start,
                          const = const,
                          location = location,
                          threshold = threshold,
                          # lme4_REML = lme4_REML,
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
                             # lme4_REML = lme4_REML,
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
                                  # lme4_REML = lme4_REML,
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
                             # lme4_REML = lme4_REML,
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
                                  # lme4_REML = lme4_REML,
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
    lambda_seq <- seq(lambda_max^(1 / lsl), 1e-3^(1 / lsl),
                      length.out = n_lambda)
    lambda_seq <- lambda_seq^lsl
  }
  
  ##############################################################################
  # calculate beta paths
  
  if(measure[1] == "CV"){
    fold_pool <- sort(table(cv_group), decreasing = TRUE)
    groups_in_folds <- lapply(
      seq_len(n_folds), function(i) as.numeric(names(fold_pool[i]))
    )
    if(length(fold_pool) > n_folds){
      for(i in seq(n_folds + 1L, length(fold_pool))){
        curr_obs <- sapply(groups_in_folds, function(vec){
          sum(cv_group %in% vec)
        })
        coi <- which.min(curr_obs)
        groups_in_folds[[coi]] <- c(
          groups_in_folds[[coi]],
          as.numeric(names(fold_pool[i]))
        )
      }
    }
    fold_i <- sapply(cv_group, function(gi){
      which(sapply(groups_in_folds, function(vec){
        gi %in% vec
      }))
    })
    
    if(method %in% c("blme", "glmmTMB")){
      mse_mat <- matrix(0, nrow = n_lambda, ncol = n_folds)
      # test_data_list <- vector(mode = "list", length = n_folds)
      # model_list <- vector(mode = "list", length = n_folds)
    } else{
      pred_mse <- numeric(n_folds)
    }
    
    for(j in seq_len(n_folds)){
      train_data <- data[fold_i != j,]
      test_data <- data[fold_i == j,]
      
      # test_data_list[[j]] <- test_data
      
      default_start <- c(
        mean(train_data[[response_name]], na.rm = TRUE),
        numeric(length(varnames))
      )
      new_start <- default_start
      
      if(method %in% c("blme", "glmmTMB")){
        for(i in seq_len(n_lambda)){
          model <- fitModel(formula = formula,
                            data = train_data,
                            method = method,
                            lambda = lambda_seq[i],
                            const = const,
                            location = location,
                            threshold = threshold,
                            max_iter = max_iter,
                            start = new_start,
                            max_attempts = max_attempts[1])
          
          if(method == "blme"){
            if(isS4(model$model)){
              pred <- predict(model$model, newdata = test_data,
                              type = "response", allow.new.levels = TRUE)
              mse_mat[i,j] <- mean(
                (test_data[[response_name]] - pred)^2,
                na.rm = TRUE
                )
            } else{
              mse_mat[i,j] <- NA
            }
          } else if(class(model$model) == "glmmTMB"){
            pred <- predict(model$model, newdata = test_data, type = "response",
                            allow.new.levels = TRUE)
            mse_mat[i,j] <- mean((test_data[[response_name]] - pred)^2,
                                 na.rm = TRUE)
          } else{
            mse_mat[i,j] <- NA
          }
          
          new_start <- svf * model$beta + (1 - svf) * default_start
          
          if(i %% message_value == 0){
            message("Fold ", j, ": ", i, " / ", n_lambda)
          }
        }
      } else{
        model <- fitModel(formula = formula,
                          data = train_data,
                          method = method,
                          const = const,
                          location = location,
                          threshold = threshold,
                          # lme4_REML = lme4_REML,
                          max_iter = max_iter,
                          start = default_start,
                          max_attempts = max_attempts[1])
        
        if(isS4(model$model)){
          pred <- predict(model$model, newdata = test_data,
                          type = "response", allow.new.levels = TRUE)
          pred_mse[j] <- mean((test_data[[response_name]] - pred)^2,
                              na.rm = TRUE)
        } else{
          pred_mse[j] <- NA
        }
      }
    }
    
    if(method %in% c("blme", "glmmTMB")){
      pred_mse <- rowMeans(mse_mat, na.rm = TRUE)
      lambda_opt <- lambda_seq[which.min(pred_mse)]
      overview <- data.frame(lambda = lambda_seq,
                             rmse = sqrt(pred_mse))
      pred_rmse <- sqrt(min(pred_mse, na.rm = TRUE))
    } else{
      pred_rmse = sqrt(mean(pred_mse, na.rm = TRUE))
      lambda_opt <- NULL
      overview <- NULL
      k <- NULL
    }
    
  } else{
    new_start <- default_start
    pred_rmse <- NA
    
    LL <- numeric(n_lambda)
    n_nonzero_param <- numeric(n_lambda)
    test_data_list <- NULL
    
    for(i in seq_len(n_lambda)){
      pen <- lambda_seq[i]
      
      model <- fitModel(formula = formula,
                        data = data,
                        method = method,
                        lambda = pen,
                        const = const,
                        location = location,
                        threshold = threshold,
                        # lme4_REML = lme4_REML,
                        max_iter = max_iter,
                        start = new_start,
                        # LL_with_RE = LL_with_RE,
                        max_attempts = max_attempts[1])
      
      new_start <- svf * model$beta + (1 - svf) * default_start
      
      LL[i] <- model$LL
      n_nonzero_param[i] <- model$nonzero
      
      if(i %% message_value == 0){
        message(i, " / ", n_lambda)
      }
    }
    
    if(measure[1] == "BIC"){
      BIC <- -2 * LL + log(nrow(data)) * n_nonzero_param
      
      lambda_opt <- lambda_seq[which.min(BIC)]
      overview <- data.frame(lambda = lambda_seq,
                             logLik = LL,
                             nonzero = n_nonzero_param,
                             BIC = BIC)
    } else{
      AIC <- -2 * LL + 2 * n_nonzero_param
      
      lambda_opt <- lambda_seq[which.min(AIC)]
      overview <- data.frame(lambda = lambda_seq,
                             logLik = LL,
                             nonzero = n_nonzero_param,
                             AIC = AIC)
    }
  }
  
  ##############################################################################
  # optimal model & output
  
  opt_model <- fitModel(formula = formula,
                        data = data,
                        method = method,
                        lambda = lambda_opt,
                        const = const,
                        location = location,
                        threshold = threshold,
                        # lme4_REML = lme4_REML,
                        max_iter = max_iter,
                        start = default_start,
                        # LL_with_RE = LL_with_RE,
                        max_attempts = max_attempts[1])
  
  return(list(
    param = opt_model$beta,
    overview = overview,
    pred_rmse = pred_rmse,
    model = opt_model,
    # test_data = test_data_list,
    lambda_max_search_runs = k
  ))
}
