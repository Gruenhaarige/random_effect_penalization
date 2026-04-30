library(lme4)
library(rpql)
library(Matrix)
library(MASS)

source("functions.R") 

{
  beta_1 <- c(5, 2)
  covX_1 <- matrix(1, nrow = 1) 
  covU_1 <- matrix(4)
  set.seed(123)
  data_1 <- generateData(
    beta = beta_1, n = 500, covX = covX_1, n_groups = 25, 
    grV = 1, slope = 0, covU = covU_1, sigma_eps = 1.5
  )
  model_1 <- lmer(y ~ x1 + (1 | G1), data = data_1$data)
  rm("beta_1", "covX_1", "covU_1")
}

y <- data_1$data$y
X <- model.matrix(~ x1, data = data_1$data)
colnames(X)[1] <- "Intercept"
Z1 <- list(model.matrix(~ 1, data = data_1$data))
colnames(Z1[[1]])[1] <- "Intercept"
id1 <- list(as.factor(data_1$data$G1))

init_weights <- function(fit_init) {
  beta_init <- fixef(fit_init)
  weights_fixed <- 1 / abs(beta_init)
  
  std_devs <- attr(VarCorr(fit_init)[[1]], "stddev")
  weights_random <- 1 / abs(std_devs)
  
  return(list(
    fixed = weights_fixed,
    random = list(weights_random) 
  ))
}

fit_init <- lmer(y ~ x1 + (1 | G1), data = data_1$data, REML = FALSE) 

(pen_weights <- init_weights(fit_init))



set.seed(123)

(lambda <- lseq(0.3, 1, length = 20))

fit <- rpqlseq(y = y, 
            X = X, 
            Z = Z1, 
            id = id1, 
            family = gaussian(), 
            pen.type = "adl",   
            pen.weights = pen_weights, 
            lambda = lambda, 
            save.data = FALSE)

summary(fit$best.fits)

fit$best.fits[[2]]$fixef
# summary(fit$best.fits[[2]]$ranef[[1]])
fit$best.fits[[2]]$ran.cov
fit$best.fits[[2]]$phi
fit$best.fits[[2]]$lambda

summary(model_1)

