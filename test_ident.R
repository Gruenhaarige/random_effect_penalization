q <- 2
ident.tilde = diag(q-1)
for (i in 2:(q-1))
{
    ident.tilde = cbind(ident.tilde, rbind(matrix(0,nrow=i-1,ncol=q-i),diag(q-i)))
}
print(ident.tilde)
