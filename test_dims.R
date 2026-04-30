n.i <- 25
q <- 2
exp.bhat <- c(0.1, 0.2)
K.matrix <- matrix(c(1, 0), nrow=2)
Z.bd <- matrix(rnorm(50), nrow=25, ncol=2)
new.lambda <- c(1.5, 0.5)
ident.tilde <- matrix(1, nrow=1, ncol=1)

A <- as.matrix(rep(1, n.i), ncol=1) %*% exp.bhat %*% K.matrix
B <- Z.bd[, 2:2] %*% diag(new.lambda[-1], nrow=1) %*% ident.tilde
print(dim(A))
print(dim(B))
A * B
