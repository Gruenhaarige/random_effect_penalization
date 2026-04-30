library(lme4)
library(MASS)
library(mvtnorm)
library(Matrix)

source("R/functions.R") 
source("R/method_lin.R")

fixed_formula <- ~ x1
random_formula <- ~ 1
grouping_var <- data_1$data$G1

data_list <- split(data_1$data, grouping_var)

y_list <- lapply(data_list, function(df) as.matrix(df$y))
x_list <- lapply(data_list, function(df) model.matrix(~ x1, data=df))
z_list <- lapply(data_list, function(df) model.matrix(~ 1, data=df))

init_fit <- lmer(y ~ x1 + (1 | G1), data = data_1$data, REML = TRUE)

vc <- VarCorr(init_fit)
sig.init <- (attributes(vc)$sc)^2
D.init <- as.matrix(vc[[1]])
beta.init <- fixef(init_fit)


lambda_re <- fit$best.fits[[2]]$lambda[2]
eps_reml <- 1e-5
lambda_fs <- fit$best.fits[[2]]$lambda[1]
eps_pco <- 1e-5

pe1 <- Pen.fs(lambda = lambda_re, 
              x = x_list, 
              y = y_list, 
              z = z_list, 
              D.init = D.init, 
              sig.init = sig.init, 
              eps = eps_reml)

final_fit <- pco(x = x_list, 
                 y = y_list, 
                 z = z_list, 
                 beta.init = pe1$beta, 
                 D.init = pe1$D, 
                 sig.init = pe1$sig, 
                 lambda = lambda_fs, 
                 eps = eps_pco)

print(final_fit$beta)
print(final_fit$D)
print(final_fit$sig)
