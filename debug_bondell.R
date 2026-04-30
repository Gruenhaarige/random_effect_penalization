source("R/generate_data.R")
source("R/method_wrappers.R")

n <- 500
p <- 8
n_signals <- 3    
n_groups <- 20
grV <- 1        
slope <- 1       
sigma_eps <- 1

covX <- diag(p)
covU <- diag(slope + 1)

beta_true <- numeric(p + 1)
signal_idx <- seq_len(n_signals)
beta_true[-1][signal_idx] <- 0.5 * (-1)^(signal_idx + 1L)
names(beta_true) <- c("(Intercept)", paste0("x", 1:p))

set.seed(123)
sim <- generateData(beta = beta_true, n = n, covX = covX, n_groups = n_groups, grV = grV, slope = slope, covU = covU, sigma_eps = sigma_eps)
data <- sim$data
data$id <- data$G1

cat("Running bondell wrapper...\n")
res <- tryCatch(run_bondell_method(data = data, random_slope_var = "x1"), error=function(e) e)
print(res)
