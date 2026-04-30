library(MASS)
library(lme4)
library(quadprog)
library(mvtnorm)
library(mlmhelpr)

source("PenLME-2cufbjn.r")

data("hsb", package = "mlmhelpr")

hsb_scaled <- hsb
for (var in c("ses", "size", "pracad", "disclim", "himinty", "meanses")) {
  hsb_scaled[[var]] <- scale(hsb_scaled[[var]])
}


fixed_formula <- ~ ses + size + pracad + disclim + himinty + meanses
random_formula <- ~ ses + minority + female + catholic

y <- as.vector(hsb_scaled$mathach)
subject <- as.numeric(as.factor(hsb_scaled$id))

X <- model.matrix(fixed_formula, data = hsb_scaled)

Z <- model.matrix(random_formula, data = hsb_scaled)
Z <- Z[, -1, drop = FALSE] 

fit_bondell <- Pen.LME(y = y, X = X, Z = Z, subject = subject, 
                       t.fracs = seq(1, 0.1, -0.1))

print(round(fit_bondell$fixed, 4))

Sigma_RE <- diag(fit_bondell$stddev) %*% fit_bondell$corr %*% diag(fit_bondell$stddev)
rownames(Sigma_RE) <- c("(Intercept)", colnames(Z))
colnames(Sigma_RE) <- c("(Intercept)", colnames(Z))
print(round(Sigma_RE, 4))

cat("\nResidual Variance (Sigma^2):\n")
print(fit_bondell$sigma.2)

cat("\nBest Tuning Parameter (t.frac):\n")
print(fit_bondell$t.frac)