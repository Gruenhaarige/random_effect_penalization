library(microbenchmark)
library(minqa)
library(ggplot2)

benchmark_results_1 <- microbenchmark(
  optim_1 = {
    result_bfgs <- optim(
      par = beta_1_initial, 
      fn = lasso_fixed_effects_objective,
      Y = Y_1,
      X_fixed = X_1,
      V_inv = V_inv_1, 
      lambda = 20,
      method = "BFGS", 
      control = list(maxit = 1000)
    )
  },
  minqa_1 = {
    result_bobyqa <- bobyqa(
      par = beta_1_initial, 
      fn = lasso_fixed_effects_objective,
      Y = Y_1,
      X_fixed = X_1,
      V_inv = V_inv_1, 
      lambda = 20,
      control = list(maxfun = 1000)
    )
  },
  times = 100
)
benchmark_results_1

benchmark_results_2 <- microbenchmark(
  optim_2 = {
    result_bfgs <- optim(
      par = beta_2_initial, 
      fn = lasso_fixed_effects_objective,
      Y = Y_2,
      X_fixed = X_2,
      V_inv = V_inv_2, 
      lambda = 20,
      method = "BFGS", 
      control = list(maxit = 2000)
    )
  },
  minqa_2 = {
    result_bobyqa <- bobyqa(
      par = beta_2_initial, 
      fn = lasso_fixed_effects_objective,
      Y = Y_2,
      X_fixed = X_2,
      V_inv = V_inv_2, 
      lambda = 20,
      control = list(maxfun = 1000)
    )
  },
  times = 100
)
benchmark_results_2


benchmark_results_2_g1 <- microbenchmark(
  optim_2_g1 = {
    result_bfgs <- optim(
      par = beta_2_g1_initial, 
      fn = lasso_fixed_effects_objective,
      Y = Y_2, 
      X_fixed = X_2, 
      V_inv = V_inv_2_g1, 
      lambda = 20,
      method = "BFGS", 
      control = list(maxit = 100)
    )
  },
  minqa_2_g1 = {
    result_bobyqa <- bobyqa(
      par = beta_2_g1_initial, 
      fn = lasso_fixed_effects_objective,
      Y = Y_2,
      X_fixed = X_2,
      V_inv = V_inv_2_g1, 
      lambda = 20,
      control = list(maxfun = 1000)
    )
  },
  times = 100
)
benchmark_results_2_g1


load("Rohfassung/benchmark_results.RData")

autoplot(benchmark_results_1) + 
  ggtitle("Benchmark: Full Model 1 (G1)") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, angle = 90))

autoplot(benchmark_results_2) + 
  ggtitle("Benchmark: Full Model 2 (G1 + G2)") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, angle = 90))

autoplot(benchmark_results_2_g1) + 
  ggtitle("Benchmark: Full Model 2 (G1)") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, angle = 90))
