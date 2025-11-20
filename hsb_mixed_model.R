library(lme4)
library(mlmhelpr)

source("R code.R")

data("hsb", package = "mlmhelpr")
hsb_scaled <- hsb

for (var in c("ses", "size", "pracad", "disclim", "himinty", "meanses")) {
  hsb_scaled[[var]] <- scale(hsb_scaled[[var]])
}

fixed_formula <- ~ ses + size + pracad + disclim + himinty + meanses
random_formula <- ~ 1 + ses + minority + female + catholic
grouping_var <- hsb_scaled$id

hsb_list <- split(hsb_scaled, grouping_var)
n_subjects <- length(hsb_list)

y_list <- lapply(hsb_list, function(df) as.matrix(df$mathach))

x_list <- lapply(hsb_list, function(df) model.matrix(fixed_formula, data=df))

z_list <- lapply(hsb_list, function(df) model.matrix(random_formula, data=df))

y1 <- do.call(rbind, y_list)
x1 <- do.call(rbind, x_list)
z1 <- do.call(rbind, z_list)
subject <- kronecker(1:n_subjects, rep(1, nrow(hsb_list[[1]]))) 


subject_ids <- rep(names(hsb_list), times = sapply(hsb_list, nrow))
subject <- as.factor(subject_ids)

lmer_data <- data.frame(
  y_response = y1,
  x_matrix = I(x1), 
  z_matrix = I(z1), 
  subject_id = subject
)

ob <- lmer(y_response ~ x_matrix - 1 + (0 + z_matrix | subject_id), data=lmer_data, REML=TRUE) 

print(summary(ob))



hh <- VarCorr(ob)
sig.init <- (attributes(hh)$sc)^2
D.init <- hh[[1]]


# STAGE 1
lambda_re <- 5
eps_reml <- 1e-5 

pe1 <- Pen.fs(lambda=lambda_re, x=x_list, y=y_list, z=z_list, 
              D.init=D.init, sig.init=sig.init, eps=eps_reml)

# pe2 <- Pen.reml(lambda=lambda_re, x=x_list, y=y_list, z=z_list,
#                D.init=pe1$D, sig.init=pe1$sig, eps=eps_reml)


# STAGE 2

lambda_fs <- 10 
eps_pco <- 1e-5

final_fit <- pco(x=x_list, y=y_list, z=z_list, 
                 beta.init=pe1$beta, 
                 D.init=pe1$D, 
                 sig.init=pe1$sig, 
                 lambda=lambda_fs, 
                 eps=eps_pco)

round(final_fit$sig, 4)
print(round(final_fit$D, 4))

beta_names <- colnames(x_list[[1]])
beta_output <- data.frame(Estimate = final_fit$beta)
rownames(beta_output) <- beta_names
print(round(beta_output, 4))


