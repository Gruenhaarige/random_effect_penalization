library(lme4)
library(mlmhelpr)
library(rpql)    
library(Matrix)

data("hsb", package = "mlmhelpr")

hsb_scaled <- hsb
for (var in c("ses", "size", "pracad", "disclim", "himinty", "meanses")) {
  hsb_scaled[[var]] <- scale(hsb_scaled[[var]])
}


fixed_formula <- ~ ses + size + pracad + disclim + himinty + meanses
random_formula <- ~ ses + minority + female + catholic

y <- hsb_scaled$mathach
X <- model.matrix(fixed_formula, data = hsb_scaled)
Z <- model.matrix(random_formula, data = hsb_scaled)

id_vec <- list(as.factor(hsb_scaled$id))

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

lmer_formula <- mathach ~ ses + size + pracad + disclim + himinty + meanses + 
  (1 + ses + minority + female + catholic | id)
fit_init <- lmer(lmer_formula, data = hsb_scaled, REML = FALSE) 

(pen_weights <- init_weights(fit_init))
names(pen_weights$fixed)[names(pen_weights$fixed) == "(Intercept)"] <- "Intercept"
names(pen_weights$random[[1]])[names(pen_weights$random[[1]]) == "(Intercept)"] <- "Intercept"
colnames(X)[colnames(X) == "(Intercept)"] <- "Intercept"
colnames(Z)[colnames(Z) == "(Intercept)"] <- "Intercept"


set.seed(123)
lambda <- lseq(1e-4, 1, length = 20) 
lambda_list <- list(lambda, lambda)
fit_hui <- rpql(y = y, 
                X = X, 
                Z = list(Z), 
                id = id_vec, 
                family = gaussian(), 
                pen.type = "adl",   
                pen.weights = pen_weights,
                lambda = lambda[1],
                save.data = FALSE) 

best_idx <- which.min(fit_hui$bic)


# Extract Coefficients for the best model
# Fixed Effects
beta_best <- fit_hui$fix.ef[best_idx, ]
cat("\nSelected Fixed Effects (Beta):\n")
print(round(beta_best[beta_best != 0], 4)) # Print non-zero

# Random Effects Covariance (D)
# rpql returns the variance components (elements of D).
# The structure is stored in $ran.cov (List of D matrices for each lambda)
D_best <- fit_hui$ran.cov[[best_idx]][[1]] # First group (id)

# Add names to D matrix for readability
colnames(D_best) <- colnames(Z_slopes)
rownames(D_best) <- colnames(Z_slopes)

cat("\nRandom Effects Covariance Matrix (D):\n")
print(round(D_best, 4))

# Residual Variance
# Stored in $phi (dispersion parameter)
cat("\nResidual Variance (Sigma^2):\n")
print(round(fit_hui$phi[best_idx], 4))

# 5. Check if 'catholic' slope was removed
# Check variance of 'catholic' in D_best
catholic_var <- D_best["catholic", "catholic"]
if (abs(catholic_var) < 1e-4) {
  cat("\nSUCCESS: 'catholic' random slope was shrunk to ZERO.\n")
} else {
  cat("\nRESULT: 'catholic' random slope was RETAINED (Var =", round(catholic_var, 4), ").\n")
}
