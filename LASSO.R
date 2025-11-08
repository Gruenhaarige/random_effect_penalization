library(MASS)
library(lme4)
library(Matrix)


## Model 1

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
  V <- sigma_u0^2 * (Z_matrix %*% t(Z_matrix)) + sigma_eps^2 * diag(n)
  V_inv <- solve(V)
  return(V_inv)
}

Y_1 <- data_1$data[, 1]
X_1 <- model.matrix(model_1)
beta_1_initial <- fixef(model_1)

V_inv_1 <- calculate_V_inv_simple(data = data_1$data, group_factor = "G1", 
                                  sigma_u0 = data.frame(VarCorr(model_1))$sdcor[1], 
                                  sigma_eps = data.frame(VarCorr(model_1))$sdcor[2])

lambda_values <- c(0, 0.1, 1, 5, 20)
results_1_lasso <- data.frame(Lambda = lambda_values, Intercept = NA, x1 = NA)

for (i in 1:length(lambda_values)) {
  lambda <- lambda_values[i]
  result <- optim(
    par = beta_m1_initial, 
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




# --------------------------------------------------------------------------------

## Model 2 

calculate_V_inv_complex <- function(data, model_fit) {
  vc <- VarCorr(model_fit)
  sigma_u_g1 <- vc$G1[1:dim(vc$G1)[1], 1:dim(vc$G1)[2]]
  sigma_u_g2 <- vc$G2[1:dim(vc$G2)[1], 1:dim(vc$G2)[2]]
  
  sigma_eps_sq <- attr(vc, "sc")^2 
  
  n <- nrow(data)
  
  Z_total <- getME(model_fit, "Z")
  
  n_g1_groups <- length(unique(data$G1))
  Sigma_u_g1_block <- bdiag(rep(list(sigma_u_g1), n_g1_groups))
  
  n_g2_groups <- length(unique(data$G2))
  Sigma_u_g2_block <- bdiag(rep(list(sigma_u_g2), n_g2_groups))
  
  Sigma_u_total <- bdiag(Sigma_u_g1_block, Sigma_u_g2_block)
  
  V <- as.matrix(Z_total %*% Sigma_u_total %*% t(Z_total) + sigma_eps_sq * Diagonal(n))
  V_inv <- solve(V)
  
  return(V_inv)
}

V_inv_2 <- calculate_V_inv_complex(data_2$data, model_2)

Y_2 <- data_2$data[, 1]
X_2 <- model.matrix(model_2)
beta_2_initial <- fixef(model_2)

lambda_values_m2 <- c(0, 0.5, 5, 20, 50)
results_2_lasso <- data.frame(Lambda = lambda_values_m2, Intercept = NA, x1 = NA, x2 = NA)

for (i in 1:length(lambda_values_m2)) {
  lambda <- lambda_values_m2[i]
  
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



## Model 2 glmmLasso

library(glmmLasso)

data_2_factor <- data_2
data_2_factor$data$G1 <- as.factor(data_2_factor$data$G1)

model_2_lasso <- glmmLasso(
  fix = y ~ x1 + x2,
  rnd = list(G1 = ~ 1 + x1, G2 = ~ 1),
  family = gaussian(link = "identity"),
  data = data_2_factor$data,
  lambda = 20  
)
model_2_lasso
results_2_lasso[6, "Lambda"] <- 20
results_2_lasso[6, 2:4] <- model_2_lasso$coefficients

results_2_lasso[c(4,6),]
