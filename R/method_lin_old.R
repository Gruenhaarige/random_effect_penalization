library("mvtnorm")
library("MASS")
library("lme4")

trace <- function(A)
{
  sum(diag(A))
}

FdD=function(di,q){
  if(q>1){
    k <- q*(q+1)/2
    dD <- matrix(0,q,q)
    DD <- rep(0,k)
    DD[di] <- 1
    dD[lower.tri(dD,diag=T)] <- DD
    dD <- dD+t(dD)-diag(diag(dD))
  }else{
    dD <- matrix(1)
  }
  list(dD)
}

Fsc <- function(r,dV,si,sii,ei,H00,Ajk){
  trace(sii%*%(ei%*%t(ei)-si)%*%sii%*%dV[[r]])+trace(H00%*%Ajk[[r]])
}

FH <- function(r,s,sii,si,ei,dV){
  trace(sii%*%dV[[r]]%*%sii%*%(2*ei%*%t(ei)-si)%*%sii%*%dV[[s]])
}

FdV <- function(r,z,dD,i){
  list(z[[i]]%*%dD[[r]]%*%t(z[[i]]))
}


FSFd2 <- function(r,s,sii,si,ei,dV){
  trace(sii%*%dV[[r]]%*%sii%*%dV[[s]])
}


Fploglik <- function(DDsig,beta,z,x,y,DDsig0,lambda,weight){
  n <- length(x)
  ni <- mapply(length,y)
  p <- ncol(x[[1]])
  q <- ncol(z[[1]])
  k <- q*(q+1)/2
  
  DD <- DDsig[1:k]
  sig <- DDsig[k+1]
  De <- matrix(0,nrow=q,ncol=q)
  De[lower.tri(De,diag=T)] <- DD
  if(ncol(De)>1){
    De <- De+t(De)-diag(diag(De))
  }
  
  ll <- 0
  rl <- 0
  for(i in 1:n){
    si <- z[[i]]%*%De%*%t(z[[i]])+sig*diag(1,ni[i])
    sii <- ginv(si)
    rl <- rl+t(x[[i]])%*%sii%*%x[[i]]
    ei <- y[[i]]-x[[i]]%*%beta
    if(det(si)<=0){return(-9e10)
    }else {des <- log(det(si))}
    
    ll <- ll-0.5*des-0.5*t(ei)%*%sii%*%ei
  }
  if(det(rl)<=0){return(-9e10)
  }else{dett <- log(det(rl))}
  
  return(ll-0.5*dett-lambda*sum(abs(DDsig*weight/DDsig0)))
  
}



FA <- function(j,x,dV,sii,i){
  list(t(x[[i]])%*%sii[[i]]%*%dV[[i]][[j]]%*%sii[[i]]%*%x[[i]])
}




FRH <- function(m,j,x,n,H00,p,dV,sii,Ajk){
  t1 <- matrix(0,p,p)
  t2 <- matrix(0,p,p)
  t3 <- matrix(0,p,p)
  for(i in 1:n){
    dercj <- dV[[i]][[j]]
    dercm <- dV[[i]][[m]]
    t1 <- t1+Ajk[[i]][[m]]
    t2 <- t2+Ajk[[i]][[j]]
    t3 <- t3-t(x[[i]])%*%sii[[i]]%*%(dercm%*%sii[[i]]%*%dercj+dercj%*%sii[[i]]%*%dercm)%*%sii[[i]]%*%x[[i]]
  }
  -trace(H00%*%t1%*%H00%*%t2)-trace(H00%*%t3)
}

Fbeta <- function(x,y,z,De,sig){
  n <- length(y)
  p <- ncol(x[[1]])
  ni <- mapply(length,y)
  bet1 <- rep(0,p)
  H11 <- matrix(0,p,p)
  for(i in 1:n){
    si <- z[[i]]%*%De%*%t(z[[i]])+sig*diag(1,ni[i])
    sii <- ginv(si)
    tss <- t(x[[i]])%*%sii
    H11 <- H11+tss%*%x[[i]]
    bet1 <- bet1+tss%*%y[[i]]
  }
  beta <- ginv(H11)%*%bet1
  beta
}

Pen.fs <- function(lambda,x,y,z,D.init,sig.init,eps){
  
  n <- length(x)
  ni <- mapply(length,y)
  p <- ncol(x[[1]])
  n.tot <- sum(ni)
  q <- ncol(z[[1]])
  q0 <- q
  k <- q*(q+1)/2
  De <- D.init
  sig <- sig.init
  beta <- Fbeta(x,y,z,De,sig)
  dD <- sapply(1:k,FdD,q)
  weight <- diag(rep(1,q))
  weight <- weight[lower.tri(weight,diag=T)]
  weight <- c(weight,0)
  DDsig0 <- c(D.init[lower.tri(D.init,diag=T)],sig.init)
  DDsignew <- DDsig0
  DD0 <- DDsig0[1:k]
  sig0 <- DDsig0[k+1]
  step <- 1
  maxstep <- 100
  record <- seq(q)
  diffll <- 10
  converge <- F
  while(step<maxstep&&diffll>1){
    DDsig <- DDsignew
    FSsc <- rep(0,1+k)
    FSH <- matrix(0,nrow=k+1,ncol=k+1)
    H0 <- matrix(0,nrow=p,ncol=p)
    
    Ajk <- list(NA)
    length(Ajk) <- n
    si <- list(NA)
    length(si) <- n
    sii <- si
    dV <- list(NA)
    length(dV) <- n
    ei <- list(NA)
    length(ei) <- n
    for(i in 1:n){
      dV[[i]] <- sapply(1:k,FdV,z,dD,i)
      dV[[i]][[k+1]] <- diag(1,ni[i])
      si[[i]] <- z[[i]]%*%De%*%t(z[[i]])+sig*diag(1,ni[i])
      sii[[i]] <- ginv(si[[i]])
      ei[[i]] <- y[[i]]-x[[i]]%*%beta
      H0 <- H0+t(x[[i]])%*%sii[[i]]%*%x[[i]]
      Ajk[[i]] <- sapply(1:(k+1),FA,x,dV,sii,i)
    }
    
    H00 <- ginv(H0)
    
    for(i in 1:n){
      FSsc <- FSsc+0.5*sapply(1:(k+1),Fsc,dV[[i]],si[[i]],sii[[i]],ei[[i]],H00,Ajk[[i]])
      H2222 <- list(NA)
      length(H2222) <- k+1
      for(j in 1:(k+1)){
        H2222[[j]] <- sapply(j:(k+1),FSFd2,j,sii[[i]],si[[i]],ei[[i]],dV[[i]])
        H2222[[j]] <- c(rep(0,(j-1)),H2222[[j]])
      }
      H222 <- do.call("cbind",H2222)
      H222 <- H222+t(H222)-diag(diag(H222))
      FSH <- FSH-0.5*H222
    }
    
    RH <- list(NA)
    length(RH) <- k+1
    for(j in 1:(k+1)){
      RH[[j]] <- sapply(j:(k+1),FRH,j,x,n,H00,p,dV,sii,Ajk)
      RH[[j]] <- c(rep(0,(j-1)),RH[[j]])
    }
    RH1 <- do.call("cbind",RH)
    RH1 <- RH1+t(RH1)-diag(diag(RH1))
    RHH <- 0.5*RH1
    
    FSsc <- FSsc-lambda*weight*sign(DDsig)/abs(DDsig0)
    FSH <- FSH-RHH-diag(lambda*weight/abs(DDsig*DDsig0))
    
    llold <- Fploglik(DDsig,beta,z,x,y,DDsig0,lambda,weight)
    llnew <- llold-1
    mm <- 1
    la <- 1
    gm <- ginv(FSH)%*%FSsc
    while(llnew<=llold&&mm<25){
      DDsignew <- DDsig-la*gm
      llnew <- Fploglik(DDsignew,beta,z,x,y,DDsig0,lambda,weight)
      la <- 1/2^mm
      mm <- mm+1
    }
    diffll <- abs(llnew-llold)
    
    signew <- DDsignew[k+1]
    DDnew <- DDsignew[1:k]
    Dnew <- matrix(0,nrow=q,ncol=q)
    Dnew[lower.tri(Dnew,diag=T)] <- DDnew
    if(ncol(Dnew)>1)Dnew <- Dnew+t(Dnew)-diag(diag(Dnew))
    
    ad <- abs(diag(Dnew))<=eps
    Dnew[ad,] <- 0
    Dnew[,ad] <- 0
    DDsignew <- c(Dnew[lower.tri(Dnew,diag=T)],signew)
    
    for(j in 1:n){
      z[[j]] <- as.matrix(z[[j]][,!ad])
    }
    record <- record[!ad]
    
    DDnew <- Dnew[lower.tri(Dnew,diag=T)]
    DD <- DDnew[abs(DDnew)>0]
    DD0 <- DD0[abs(DDnew)>0]
    De <- Dnew[!ad,!ad]
    De <- as.matrix(De)
    sig <- signew
    beta <- Fbeta(x,y,z,De,sig)
    q <- dim(De)[2]
    k <- q*(q+1)/2
    DDsig0 <- c(DD0,sig0)
    
    dD <- sapply(1:k,FdD,q)
    weight <- diag(rep(1,q))
    weight <- weight[lower.tri(weight,diag=T)]
    weight <- c(weight,0)
    
    DDsignew <- c(De[lower.tri(De,diag=T)],signew)
    
    step <- step+1
  }
  Df <- matrix(0,q0,q0)
  Df[record,record]=De
  fit <- NULL
  fit$beta <- beta
  fit$D <- Df
  fit$sig <- sig
  return(fit)
}



Pen.reml <- function(lambda,x,y,z,D.init,sig.init,eps){
  
  n <- length(x)
  ni <- mapply(length,y)
  p <- ncol(x[[1]])
  n.tot <- sum(ni)
  q <- ncol(z[[1]])
  q0 <- q
  k <- q*(q+1)/2
  De <- D.init
  sig <- sig.init
  beta <- Fbeta(x,y,z,De,sig)
  dD <- sapply(1:k,FdD,q)
  weight <- diag(rep(1,q))
  weight <- weight[lower.tri(weight,diag=T)]
  weight <- c(weight,0)
  DDsig0 <- c(D.init[lower.tri(D.init,diag=T)],sig.init)
  DDsignew <- DDsig0
  DD0 <- DDsig0[1:k]
  sig0 <- DDsig0[k+1]
  
  step <- 1
  maxstep <- 100
  record <- seq(q)
  converge <- F
  while(converge==F&&step<maxstep&&ncol(De)>1){
    
    DDsig <- DDsignew
    H0 <- matrix(0,nrow=p,ncol=p)
    sc <- rep(0,1+k)
    H <- matrix(0,nrow=k+1,ncol=k+1)
    
    Ajk <- list(NA)
    length(Ajk) <- n
    si <- list(NA)
    length(si) <- n
    sii <- si
    dV <- list(NA)
    length(dV) <- n
    ei <- list(NA)
    length(ei) <- n
    for(i in 1:n){
      dV[[i]] <- sapply(1:k,FdV,z,dD,i)
      dV[[i]][[k+1]] <- diag(1,ni[i])
      si[[i]] <- z[[i]]%*%De%*%t(z[[i]])+sig*diag(1,ni[i])
      sii[[i]] <- ginv(si[[i]])
      ei[[i]] <- y[[i]]-x[[i]]%*%beta
      H0 <- H0+t(x[[i]])%*%sii[[i]]%*%x[[i]]
      Ajk[[i]] <- sapply(1:(k+1),FA,x,dV,sii,i)
    }
    
    H00 <- ginv(H0)
    
    for(i in 1:n){
      scb <- sapply(1:(k+1),Fsc,dV[[i]],si[[i]],sii[[i]],ei[[i]],H00,Ajk[[i]])
      sc <- sc+0.5*scb
      H2222 <- list(NA)
      length(H2222) <- k+1
      for(j in 1:(k+1)){
        H2222[[j]] <- sapply(j:(k+1),FH,j,sii[[i]],si[[i]],ei[[i]],dV[[i]])
        H2222[[j]] <- c(rep(0,(j-1)),H2222[[j]])
      }
      H222 <- do.call("cbind",H2222)
      H222 <- H222+t(H222)-diag(diag(H222))
      H <- H-0.5*H222
    }
    
    
    RH <- list(NA)
    length(RH) <- k+1
    for(j in 1:(k+1)){
      RH[[j]] <- sapply(j:(k+1),FRH,j,x,n,H00,p,dV,sii,Ajk)
      RH[[j]] <- c(rep(0,(j-1)),RH[[j]])
    }
    RH1 <- do.call("cbind",RH)
    RH1 <- RH1+t(RH1)-diag(diag(RH1))
    RHH <- 0.5*RH1
    
    sc <- sc-lambda*weight*sign(DDsig)/abs(DDsig0)
    H <- H-RHH-diag(lambda*weight/abs(DDsig*DDsig0))
    
    llold <- Fploglik(DDsig,beta,z,x,y,DDsig0,lambda,weight)
    llnew <- llold-1
    mm <- 1
    la <- 1
    
    gH <- ginv(H)%*%sc
    while(llnew<=llold&&mm<15){
      DDsignew <- DDsig-la*gH
      llnew <- Fploglik(DDsignew,beta,z,x,y,DDsig0,lambda,weight)
      la <- 1/2^mm
      mm <- mm+1
    }
    
    signew <- DDsignew[k+1]
    DDnew <- DDsignew[1:k]
    Dnew <- matrix(0,nrow=q,ncol=q)
    Dnew[lower.tri(Dnew,diag=T)] <- DDnew
    
    if(ncol(Dnew)>1)Dnew <- Dnew+t(Dnew)-diag(diag(Dnew))
    
    ad <- abs(diag(Dnew))<=eps
    Dnew[ad,] <- 0
    Dnew[,ad] <- 0
    DDsignew <- c(Dnew[lower.tri(Dnew,diag=T)],signew)
    
    for(j in 1:n){
      z[[j]] <- as.matrix(z[[j]][,!ad])
    }
    record <- record[!ad]
    if(sum((DDsignew-DDsig)^2)<eps)converge <- TRUE
    DDnew <- Dnew[lower.tri(Dnew,diag=T)]
    DD <- DDnew[abs(DDnew)>0]
    DD0 <- DD0[abs(DDnew)>0]
    De <- Dnew[!ad,!ad]
    De <- as.matrix(De)
    DDsig0 <- c(DD0,sig0)
    
    sig <- signew
    beta <- Fbeta(x,y,z,De,sig)
    q <- dim(De)[2]
    k <- q*(q+1)/2
    dD <- sapply(1:k,FdD,q)
    weight <- diag(rep(1,q))
    weight <- weight[lower.tri(weight,diag=T)]
    weight <- c(weight,0)
    DDsignew <- c(De[lower.tri(De,diag=T)],sig)
    
    step <- step+1
  }
  
  Df <- matrix(0,q0,q0)
  Df[record,record]=De
  fit <- NULL
  fit$beta <- beta
  fit$D <- Df
  fit$sig <- sig
  return(fit)
}



pco <- function(x,y,z,beta.init,D.init,sig.init,lambda,eps){
  
  D <- D.init
  sig <- sig.init
  beta <- beta.init
  
  p <- ncol(x[[1]])
  n <- length(y)
  ni <- mapply(length,y)
  
  si <- list(NA)
  length(si) <- n
  sii <- si
  siix <- si
  siixy <- matrix(0,nrow=p,ncol=1)
  siixx <- matrix(0,nrow=p,ncol=p)
  for(i in 1:n){
    si[[i]] <- z[[i]]%*%D%*%t(z[[i]])+sig*diag(1,ni[i])
    sii[[i]] <- ginv(si[[i]])
    siix[[i]] <- t(x[[i]])%*%sii[[i]]
    siixy <- siix[[i]]%*%y[[i]]+siixy
    siixx <- siix[[i]]%*% x[[i]]+siixx
  }
  
  err <- 100
  step <- 1
  maxstep <- 100
  while(err>eps&&step<maxstep){
    beta.old <- beta
    for(j in 1:p){
      S0 <- siixy[j]-siixx[j,-j]%*%beta[-j]
      if(S0>0&&S0>lambda/abs(beta.init[j])){beta[j] <- (S0-lambda/abs(beta.init[j]))/siixx[j,j]
      }else if(S0<0&&lambda/abs(beta.init[j])<abs(S0)){beta[j] <- (S0+lambda/abs(beta.init[j]))/siixx[j,j]
      }else{
        beta[j] <- 0
      }
    }
    err <- sum((beta.old-beta)^2)
    step <- step+1
  }
  fit <- NULL
  fit$beta <- beta
  fit$D <- D
  fit$sig <- sig
  return(fit)
}


