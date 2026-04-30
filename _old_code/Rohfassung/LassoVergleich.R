model_2_g1 <- lmer(y ~ x1 + x2 + (1 + x1 | G1), data = data_2$data)

calculate_V_inv_general <- function(data, model_fit) {
  vc <- VarCorr(model_fit)
  sigma_eps_sq <- attr(vc, "sc")^2 
  n <- nrow(data)
  Z_total <- getME(model_fit, "Z")
  
  sigma_u_g1 <- vc$G1[1:dim(vc$G1)[1], 1:dim(vc$G1)[2]]
  Sigma_u_total <- bdiag(bdiag(rep(list(sigma_u_g1), n_groups_i)))
  
  V <- as.matrix(Z_total %*% Sigma_u_total %*% t(Z_total) + sigma_eps_sq * Diagonal(n))
  V_inv <- solve(V)
  return(V_inv)
}

V_inv_2_g1 <- calculate_V_inv_general(data = data_2$data, model_fit = model_2_g1)
beta_2_g1_initial <- fixef(model_2_g1)

result_manual_2_g1 <- optim(
  par = beta_2_g1_initial, 
  fn = lasso_fixed_effects_objective,
  Y = Y_2, X_fixed = X_2, V_inv = V_inv_2_g1, lambda = 20,
  method = "BFGS", control = list(maxit = 2000)
)
print(result_manual_2_g1$par)




library(glmmLasso)

data_2_factor <- data_2
data_2_factor$data$G1 <- as.factor(data_2_factor$data$G1)

fit_glmmLasso_2_g1 <- glmmLasso(
  fix = y ~ x1 + x2,     
  rnd = list(G1 = ~ 1 + x1),
  data = data_2_factor$data, 
  lambda = 10,
  family = gaussian(link = "identity"), 
  control = list(print.iter = FALSE, crit.max = 100)
)

results_2_g1_lasso <- data.frame(type = c("Manuell", "With glmmLasso"), Intercept = NA, x1 = NA, x2 = NA)
results_2_g1_lasso[1, 2:4] <- result_manual_2_g1$par
results_2_g1_lasso[2, 2:4] <- model_2_lasso$coefficients

results_2_g1_lasso
