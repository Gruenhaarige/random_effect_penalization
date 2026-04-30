fitModel <- function(formula,
                     data,
                     method,
                     lambda = NULL,
                     const = 1e-5,
                     location = 0,
                     threshold = 1e-3,
                     common_scale = TRUE,
                     max_iter = 1e5,
                     start = NULL,
                     vary_lambda = FALSE,
                     max_attempts = 10L){
  
  ##############################################################################
  
  errorFun <- function(e){
    message(e)
    return(NA)
  }
  
  response_name <- all.vars(formula[[2]])
  re_obj <- lme4::findbars(formula)
  group_var <- unique(vapply(re_obj, function(x) deparse(x[[3]]), character(1)))
  
  covar <- attr(terms(nobars(formula)), "term.labels")
  param_length <- length(covar) + 1L
  param_names <- c("(Intercept)", covar)
  
  data_means <- sapply(data[covar], mean, na.rm = TRUE)
  data_sd <- sapply(data[covar], sd, na.rm = TRUE)
  
  data[covar] <- scale(data[covar])
  
  ##############################################################################
  
  if(is.null(start)){
    response_name <- all.vars(formula)[1]
    start <- c(mean(data[[response_name]]), numeric(param_length))
  }
  if(is.null(names(start))){
    names(start) <- c("(Intercept)", covar)
  }
  
  attempts <- 0L
  
  repeat{
    if(attempts == max_attempts){
      fit <- NA
      param <- NA
      LL <- NA
      n_param_total <- NA
      message("fit failed ", max_attempts, " times!")
      break
    }
    
    if(method == "glmmTMB"){
      prior_string <- paste0("lasso(",
                             location, ", ",
                             lambda, ", ",
                             const, ", ",
                             as.integer(common_scale), ")")
      prior <- data.frame(prior = prior_string,
                          class = c("fixef"),
                          coef = c(""))
      startlist = list(beta = start)
      fit <- tryCatch(glmmTMB(formula,
                              data = data,
                              priors = prior,
                              start = startlist,
                              control = glmmTMBControl(
                                rank_check = "skip",
                                optCtrl = list(eval.max = max_iter,
                                               iter.max = max_iter),
                                profile=TRUE
                              )
      ), error = function(e) errorFun(e))
      
      if(length(fit) == 1 && is.na(fit)){
        start <- start + rnorm(param_length)
        if(vary_lambda){
          lambda <- lambda + rexp(1L, rate = 20)
        }
        attempts <- attempts + 1L
        message("model fit failed ", attempts, " times!")
      } else{
        param <- fit$fit$par[names(fit$fit$par) == "beta"]
        param <- ifelse(abs(param) < threshold, 0, param)
        n_param_total <- glmmTMB:::npar(fit)
        
        # hack to exclude the log prior from the logLik output
        # see https://github.com/glmmTMB/glmmTMB/issues/1116
        pars <- fit$obj$env$parList()
        pars2 <- pars[lengths(pars)>0 & names(pars)!="b"]
        mapvec <- lapply(pars2, \(x) factor(rep(NA_real_, length(x))))
        
        fit_for_LL <- tryCatch(update(
          fit, priors = NULL, start = pars2, map = mapvec,
          control = glmmTMBControl(rank_check = "skip", profile = FALSE)
        ), error = function(e) NA)
        
        if(length(fit_for_LL) == 1 && is.na(fit_for_LL)){
          LL <- NA
        } else{
          LL <- logLik(fit_for_LL)
        }
        break
      }
      
    } else if(method == "blme"){
      fit <- tryCatch(
        blmer(formula,
              data = data,
              REML = FALSE,
              fixef.prior = lasso(mean = 0,
                                  lambda = lambda,
                                  c = const,
                                  common.scale = common_scale),
              # start = list(fixef = start), # buggy
              control = lmerControl(optimizer = "bobyqa",
                                    optCtrl = list(maxfun = max_iter))
        ), error = function(e) errorFun(e))
      
      if(!isS4(fit) && length(fit) == 1 && is.na(fit)){
        start <- start + rnorm(param_length)
        if(vary_lambda){
          lambda <- lambda + rexp(1L, rate = 20)
        }
        attempts <- attempts + 1L
        message("model fit failed ", attempts, " times!")
      } else{
        param <- fit@beta
        n_param_total <- lme4:::npar.merMod(fit)
        LL <- logLik(fit)
        param <- ifelse(abs(param) < threshold, 0, param)
        break
      }
    } else if(method == "lme4"){
      fit <- tryCatch(lmer(formula,
                           data = data,
                           start = list(fixef = start)),
                      error = function(e) errorFun(e))
      
      if(!isS4(fit) && length(fit) == 1 && is.na(fit)){
        start <- start + rnorm(param_length)
        attempts <- attempts + 1L
        message("model fit failed ", attempts, " times!")
      } else{
        param <- fit@beta
        n_param_total <- lme4:::npar.merMod(fit)
        LL <- logLik(fit)
        break
      }
    } else{
      stop("Unknown method!")
    }
  }
  
  if(!(length(param) == 1 && is.na(param))){
    names(param) <- param_names
  }
  
  return(list(
    beta = param,
    model = fit,
    LL = LL,
    # parameters shrunken to zero don't count
    nonzero = n_param_total - sum(param == 0),
    lambda = lambda,
    scale_val = list(means = data_means, sds = data_sd)
  ))
}
