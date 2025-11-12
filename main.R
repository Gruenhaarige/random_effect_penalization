library(MASS)
library(lme4)
library(glmmLasso)

source("functions.R")

## Daten und Modelle erstellen 

# Modell 1
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
  print(summary(model_1))
  rm("beta_1", "covX_1", "covU_1")
}

# Modell 2 (mit und ohne G2)
{
  beta_2 <- c(10, 3, -1)
  covX_2 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2) 
  covU_2 <- matrix(c(9, 1, 1, 4), nrow = 2) 
  set.seed(123)
  data_2 <- generateData(
    beta = beta_2, n = 1000, covX = covX_2, n_groups = 40, 
    grV = 2, slope = 1, covU = covU_2, sigma_eps = 2.0
  )
  model_2 <- lmer(y ~ x1 + x2 + (1 + x1 | G1) + (1 | G2), data = data_2$data)
  model_2_g1_only <- lmer(y ~ x1 + x2 + (1 + x1 | G1), data = data_2$data)
  print(summary(model_2))
  print(summary(model_2_g1_only))
  rm("beta_2", "covU_2", "covX_2")
}


## Lasso für die festen Effecte

# Modell 1
{
  Y_1 <- data_1$data[, 1]
  X_1 <- model.matrix(model_1)
  beta_1_initial <- fixef(model_1)
  
  V_inv_1 <- calculate_V_inv_simple(data = data_1$data, group_factor = "G1", 
                                    sigma_u0 = data.frame(VarCorr(model_1))$sdcor[1], 
                                    sigma_eps = data.frame(VarCorr(model_1))$sdcor[2])
}

{
  lambda_values <- c(0, 1, 50, 100, 500)
  results_1_lasso <- data.frame(Lambda = lambda_values, Intercept = NA, x1 = NA)
  
  for (i in 1:length(lambda_values)) {
    lambda <- lambda_values[i]
    result <- optim(
      par = beta_1_initial, 
      fn = lasso_fixed_effects_objective,
      Y = Y_1,
      X_fixed = X_1,
      V_inv = V_inv_1, 
      lambda = lambda,
      method = "BFGS", 
      control = list(maxit = 1000)
    )
    results_1_lasso[i, "Intercept"] <- result$par[1]
    results_1_lasso[i, "x1"] <- result$par[2]
  }
  print(round(results_1_lasso, 4))
  # rm("Y_1", "X_1", "beta_1_initial", "V_inv_1", "lambda_values", "i", "lambda", "result")
}

# Modell 2
{
  Y_2 <- data_2$data[, 1]
  X_2 <- model.matrix(model_2)
  beta_2_initial <- fixef(model_2)
  V_inv_2 <- calculate_V_inv_complex(data = data_2$data, model_fit = model_2, g2 = TRUE)
}

{
  lambda_values_2 <- c(0, 0.5, 5, 20, 50)
  results_2_lasso <- data.frame(Lambda = lambda_values_2, Intercept = NA, x1 = NA, x2 = NA)
  
  for (i in 1:length(lambda_values_2)) {
    lambda <- lambda_values_2[i]
    result <- optim(
      par = beta_2_initial, 
      fn = lasso_fixed_effects_objective,
      Y = Y_2,
      X_fixed = X_2,
      V_inv = V_inv_2, 
      lambda = lambda,
      method = "BFGS", 
      control = list(maxit = 2000)
    )
    results_2_lasso[i, 2:4] <- result$par
  }
  print(round(results_2_lasso, 4))
}

## Modell 2 (ohne G2)

{
  V_inv_2_g1 <- calculate_V_inv_complex(data = data_2$data, model_fit = model_2_g1_only, g2 = FALSE)
  beta_2_g1_initial <- fixef(model_2_g1_only)
}

{
  result_manual_2_g1 <- optim(
    par = beta_2_g1_initial, 
    fn = lasso_fixed_effects_objective,
    Y = Y_2, X_fixed = X_2, V_inv = V_inv_2_g1, lambda = 20,
    method = "BFGS", control = list(maxit = 100)
  )
  
  data_2_factor <- data_2$data
  data_2_factor$G1 <- as.factor(data_2_factor$G1)
  
  fit_glmmLasso_2_g1 <- glmmLasso(
    fix = y ~ x1 + x2,     
    rnd = list(G1 = ~ 1 + x1),
    data = data_2_factor, 
    lambda = 20,
    family = gaussian(link = "identity"), 
    control = list(print.iter = FALSE, crit.max = 100)
  )
  
  results_2_g1_lasso <- data.frame(type = c("Manuell", "glmmLasso"), 
                                   Intercept = NA, x1 = NA, x2 = NA)
  results_2_g1_lasso[1, 2:4] <- result_manual_2_g1$par
  results_2_g1_lasso[2, 2:4] <- fit_glmmLasso_2_g1$coefficients
  
  results_2_g1_lasso
}

