library(MuMIn)

r2 <- function(fit_object, x_list, z_list, y_list) {
  
  beta <- fit_object$beta
  D <- fit_object$D
  sigma2_eps <- fit_object$sig
  X_full <- do.call(rbind, x_list)
  pred_fixed <- X_full %*% beta
  sigma2_f <- var(as.vector(pred_fixed))
  var_random <- sapply(z_list, function(Z_i) {
    mean(diag(Z_i %*% D %*% t(Z_i)))
  })
  sigma2_r <- mean(var_random)
  
  total_var <- sigma2_f + sigma2_r + sigma2_eps
  
  r2_marginal <- sigma2_f / total_var
  r2_conditional <- (sigma2_f + sigma2_r) / total_var
  
  return(c(R2_Marginal = r2_marginal, R2_Conditional = r2_conditional, 
           Var_Fixed = sigma2_f, Var_Random = sigma2_r, Var_Resid = sigma2_eps))
}

r2_full <- r2(final_fit, x_list, z_list, y_list)
print(round(r2_full, 4))

r2_stage2 <- r2(final_fit_stage2, x_list, z_list, y_list)
print(round(r2_stage2, 4))

r.squaredGLMM(ob)

comparison_df <- data.frame(
  Model = c("Unpenalized", "Full_Penalized", "Fixed_Effects_Penalized"),
  R2m_Fixed = c(r.squaredGLMM(ob)[1], r2_full["R2_Marginal"], r2_stage2["R2_Marginal"]),
  R2c_Total = c(r.squaredGLMM(ob)[2], r2_full["R2_Conditional"], r2_stage2["R2_Conditional"])
)

comparison_df
