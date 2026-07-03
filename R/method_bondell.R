
##################################################################
#
#
#   Function to fit Penalized Mixed Model method of 
#        Bondell, Krishna, and Ghosh (2009)
#
#
#
#####################################################################



Pen.LME.fit = function(y, X, Z, subject, t.fracs, eps = 10^(-4)) {
	# Assumes Z already has intercept prepended!
	n.i = tabulate(subject)
	n.tot = sum(n.i)
	n = length(n.i)
	p = ncol(X)
	q = ncol(Z)
	
	eps.tol = 0.00000001
	
	# Initial fit to get starting values
	init.fit = tryCatch({
		lmer(y ~ X - 1 + (Z - 1 | subject))
	}, error = function(e) {
		# fallback control parameters if default optimizer fails
		lmer(y ~ X - 1 + (Z - 1 | subject), control = lmerControl(optimizer = "Nelder_Mead"))
	})
	
	est = VarCorr(init.fit)
	sigma.hat = (attributes(est)$sc)^2
	beta.hat = as.matrix(fixef(init.fit))
	beta.hatp = abs(beta.hat)
	beta.p = t(1/rbind(beta.hatp,beta.hatp))
	D.lme = as.matrix(est$subject)/sigma.hat
	
	junk = t(chol(D.lme+eps.tol*diag(q)))
	lambda.hat = diag(junk)
	lambda.hat = pmax(lambda.hat, eps.tol)
	gamma.init = diag(as.vector(1/lambda.hat))%*%junk
	gamma.hat = gamma.init[lower.tri(gamma.init)]
	lambda.p = t(as.matrix(1/lambda.hat))
	
	Z.bd = matrix(0,nrow=(n.tot),ncol=(n*q))
	W.bd = matrix(0,nrow=(n*q),ncol=(n*q))
	start.point = 1
	for (i in 1:n)
	{
		end.point = start.point + (n.i[i] - 1) 
		Z.bd[start.point:end.point,(q*(i-1)+1):(q*i)] = Z[start.point:end.point,]
		start.point = end.point + 1
	}
	W.bd = t(Z.bd)%*%Z.bd
	
	X.star = cbind(X,-X)
	X.star.quad = t(X.star)%*%X.star
	
	A.trans = rbind(diag(2*p + q), -c(beta.p,lambda.p))
	cr.full.k = kronecker(rep(1,n),diag(q))
	
	ident.tilde = diag(q-1)
	K.matrix = NULL
	for (i in 1:(q-1))
	{
		for (j in 1:(q-i))
		{
			K.matrix = cbind(K.matrix, c(rep(0,i-1),1,rep(0,q-i-1)))
		}
	}
	K.matrix=rbind(K.matrix,rep(0,q*(q-1)/2))
	if (q > 2) {
		for (i in 2:(q-1))
		{
			ident.tilde = cbind(ident.tilde, rbind(matrix(0,nrow=i-1,ncol=q-i),diag(q-i)))
		}
	}
	
	# Warm start states
	curr.beta = beta.hat
	curr.lambda = lambda.hat
	curr.gamma = gamma.hat
	curr.sigma.2 = as.numeric(sigma.hat)
	
	results = list()
	
	for (k in seq_along(t.fracs))
	{
		frac = t.fracs[k]
		new.beta = curr.beta
		new.lambda = curr.lambda
		new.gamma = curr.gamma
		sigma.2.current = curr.sigma.2
		
		t.bound = frac*(p+q)
		b.0 = c(rep(0,2*p + q), -t.bound)
		outer.converge = FALSE
		n.iter = 0
		
		while ((outer.converge==FALSE) && (n.iter < 200))
		{
			n.iter = n.iter + 1
			beta.current = beta.iterate = new.beta
			lambda.current = new.lambda
			gamma.current = new.gamma
			gamma.mat.current = diag(q)
			gamma.mat.current[lower.tri(gamma.mat.current)] = gamma.current
			full.gamma.mat = kronecker(diag(n),gamma.mat.current)
			
			n.iter1 = 0
			inner.converge = FALSE
			while ((inner.converge==FALSE) && (n.iter1 < 100))
			{	
				beta.current = new.beta
				lambda.current = new.lambda
				n.iter1 = n.iter1 + 1
				resid.vec.current = y-(X%*%beta.current)
				full.gamma.mat = kronecker(diag(n),gamma.mat.current)
				
				full.D.mat = kronecker(diag(n),diag(as.vector(lambda.current)))
				Cov.mat.temp = as.matrix(Z.bd%*%full.D.mat%*%full.gamma.mat)
				sigma.2.current = as.numeric(t(resid.vec.current)%*%solve(Cov.mat.temp%*%t(Cov.mat.temp)+diag(n.tot))%*%resid.vec.current/n.tot)
				
				full.inv.Cov.mat = solve(t(Cov.mat.temp)%*%Cov.mat.temp + diag(n*q))
				exp.bhat = full.inv.Cov.mat%*%t(Cov.mat.temp)%*%resid.vec.current
				exp.Uhat = full.inv.Cov.mat*sigma.2.current
				exp.Ghat = exp.Uhat + exp.bhat%*%t(exp.bhat)
				
				right.side.mat = as.matrix(Z.bd%*%diag(as.vector(full.gamma.mat%*%exp.bhat))%*%cr.full.k)
				lower.diag.mat = as.matrix(t(cr.full.k)%*%(W.bd * (full.gamma.mat%*%exp.Ghat%*%t(full.gamma.mat)))%*%cr.full.k)
				
				full.right.side = as.matrix(t(X.star)%*%right.side.mat)
				D.quadratic.prog = rbind(cbind(as.matrix(X.star.quad), full.right.side),cbind(t(full.right.side), lower.diag.mat))
				d.linear.prog = as.vector(t(y)%*%cbind(as.matrix(X.star), right.side.mat))
				D.quadratic.prog = D.quadratic.prog+eps.tol*diag(nrow(D.quadratic.prog))
				
				beta.lambda = solve.QP(D.quadratic.prog, d.linear.prog, t(A.trans), bvec=b.0)
				new.beta = round(beta.lambda$solution[1:p]-beta.lambda$solution[(p+1):(2*p)],6)
				new.lambda = round(beta.lambda$solution[-(1:(2*p))],6)
				
				diff = abs(beta.current-new.beta)
				if (max(c(diff))<eps) 
				{
					inner.converge = TRUE
				}
			}
			
			E.A = NULL
			start.point = 1
			d.d.t = new.lambda%*%t(new.lambda)
			
			full.A.t.A.matrix = 0*diag(q*(q-1)/2)
			T.vec = rep(0,q*(q-1)/2)
			
			for (i in 1:n)
			{
				end.point = start.point + (n.i[i] - 1) 
				E.Ai = as.matrix((as.matrix(rep(1,n.i[i]),ncol=1)%*%exp.bhat[(q*(i-1)+1):(q*i)]%*%K.matrix)*(Z.bd[start.point:end.point,(q*(i-1)+2):(q*i)]%*%diag(new.lambda[-1], nrow=length(new.lambda[-1]))%*%ident.tilde))
				E.A = rbind(E.A, E.Ai)
				start.point = end.point + 1
				
				G.i = exp.Ghat[(q*(i-1)+1):(q*i),(q*(i-1)+1):(q*i)]
				Z.Z.i = W.bd[(q*(i-1)+1):(q*i),(q*(i-1)+1):(q*i)]
				B.matrix = as.matrix(Z.Z.i * d.d.t)
				A.t.A.matrix = NULL
				Cross.matrix = NULL
				
				for (j in 1:(q-1))
				{
					U.j.matrix = diag(q)[-(1:j),]
					Cross.matrix = rbind(Cross.matrix, U.j.matrix%*%B.matrix)
					A.t.A.row = NULL
					for (k2 in 1:(q-1))
					{
						V.k.matrix = diag(q)[,-(1:k2)]
						A.t.A.row = cbind(A.t.A.row, U.j.matrix%*%B.matrix%*%V.k.matrix*G.i[j,k2])  
					}
					A.t.A.matrix = rbind(A.t.A.matrix, A.t.A.row)
				}
				T.vec = T.vec + as.vector(((t(K.matrix)%*%G.i)*Cross.matrix)%*%rep(1,q))
				full.A.t.A.matrix = full.A.t.A.matrix + A.t.A.matrix
			}
			
			A.eigen = eigen(full.A.t.A.matrix)
			A.eigen.vals = round(A.eigen$values,5)
			A.eigen.vecs = A.eigen$vectors
			eig.A = A.eigen.vals^(-1)
			eig.A[is.infinite(eig.A)] = 0
			
			A.t.A.inv = (A.eigen.vecs%*%diag(eig.A, nrow=length(eig.A))%*%t(A.eigen.vecs))
			lin.term = t(E.A)%*%resid.vec.current-T.vec
			new.gamma = round(A.t.A.inv%*%lin.term,6)
			counter.1 = 1
			for (j in 1:(q-1))
			{
				for (k2 in (j+1):q)
				{
					new.gamma[counter.1] = new.gamma[counter.1]*(new.lambda[j]>0)
					counter.1 = counter.1 + 1
				}
			}
			
			diff = abs(beta.iterate-new.beta)
			if (max(c(diff))<eps) 
			{
				outer.converge = TRUE
			}
		}
		
		# Save this fraction's state
		results[[k]] = list(
			beta = new.beta,
			lambda = new.lambda,
			gamma = new.gamma,
			sigma.2 = sigma.2.current
		)
		
		# Update warm start state for next fraction
		curr.beta = new.beta
		curr.lambda = new.lambda
		curr.gamma = new.gamma
		curr.sigma.2 = sigma.2.current
	}
	
	return(results)
}

Pen.LME = function(y, X, Z, subject, t.fracs = seq(1,0.05,-0.05), eps = 10^(-4))
{
	require(MASS)
	require(lme4)
	require(quadprog)
	require(mvtnorm)
	t.fracs = sort(t.fracs,decreasing=T)
	if (min(t.fracs) <=0) {return(cat("ERROR: All values for t.fracs must be > 0. \n"))}
	if (max(t.fracs) >1) {return(cat("ERROR: All values for t.fracs must be < 1. \n"))}
	if (eps <=0) {return(cat("ERROR: Eps must be > 0. \n"))}
	
	n.i = tabulate(subject)
	n.tot = sum(n.i)
	n.subj = length(n.i)
	Z = as.matrix(Z,nrow=n.tot)
	Z = cbind(rep(1,n.tot), Z)
	p = ncol(X)
	q = ncol(Z)	
	
	if (qr(X)$rank < p) {return(cat("ERROR: Design matrix for fixed effects is not full rank. \n"))}
	if (qr(Z)$rank < q) {return(cat("ERROR: Design matrix for random effects is not full rank. \n"))}

	# Single warm-started sweep across the t.frac grid (largest -> smallest,
	# see Pen.LME.fit's curr.* state). No cross-validation: Pen.LME is a single
	# joint model over fixed AND random effects, so BIC is computed once
	# against the full joint marginal likelihood below, using total N
	# throughout (fixed effects, random effects, and residual variance all
	# come from the same model fit on all n.tot observations).
	full_fits <- Pen.LME.fit(y, X, Z, subject, t.fracs, eps)

	# Block-diagonal Z, needed both for BIC evaluation and BLUP extraction.
	Z.bd = matrix(0,nrow=(n.tot),ncol=(n.subj*q))
	start.point = 1
	for (i in 1:n.subj)
	{
		end.point = start.point + (n.i[i] - 1)
		Z.bd[start.point:end.point,(q*(i-1)+1):(q*i)] = Z[start.point:end.point,]
		start.point = end.point + 1
	}

	# BIC selection across the grid of fractions, using the true joint
	# marginal log-likelihood (dmvnorm over y with the full Z D Z' + sigma^2*I
	# covariance), matching method_bondell_old.R's formula -- NOT the
	# independent-Gaussian approximation this replaces, which ignored the
	# random-effect covariance entirely.
	BIC.values = sapply(full_fits, function(fit) {
		beta = fit$beta
		lambda = fit$lambda
		gamma = fit$gamma
		sigma2 = fit$sigma.2

		gamma.mat = diag(q)
		gamma.mat[lower.tri(gamma.mat)] = gamma

		full.gamma.mat = kronecker(diag(n.subj), gamma.mat)
		full.D.mat = kronecker(diag(n.subj), diag(as.vector(lambda)))
		Cov.mat.temp = as.matrix(Z.bd %*% full.D.mat %*% full.gamma.mat)
		Full.cov.mat = Cov.mat.temp %*% t(Cov.mat.temp) + diag(n.tot)
		Full.Cov.est = sigma2 * Full.cov.mat
		Mean.est = X %*% beta

		loglik = -2 * dmvnorm(as.vector(y), Mean.est, Full.Cov.est, log = TRUE)
		df.par = sum(beta != 0) + sum(lambda != 0) * ((sum(lambda != 0) + 1) / 2)
		loglik + df.par * log(n.tot)
	})
	min.BIC.idx = which.min(BIC.values)
	opt.frac = t.fracs[min.BIC.idx]
	# Extract optimal fit
	opt_fit = full_fits[[min.BIC.idx]]
	beta.opt = opt_fit$beta
	lambda.opt = opt_fit$lambda
	gamma.opt = opt_fit$gamma
	sigma.2.opt = opt_fit$sigma.2

	gamma.BIC.mat = diag(q)
	gamma.BIC.mat[lower.tri(gamma.BIC.mat)] = gamma.opt
	temp.mat = diag(as.vector(lambda.opt))%*%gamma.BIC.mat
	Cov.Mat.RE = sigma.2.opt * temp.mat %*% t(temp.mat)

	# Compute optimal random effects BLUPs (exp.bhat) on the full dataset
	full.gamma.mat = kronecker(diag(n.subj), gamma.BIC.mat)
	full.D.mat = kronecker(diag(n.subj), diag(as.vector(lambda.opt)))
	Cov.mat.temp = as.matrix(Z.bd %*% full.D.mat %*% full.gamma.mat)
	resid.vec.current = y - (X %*% beta.opt)
	full.inv.Cov.mat = solve(t(Cov.mat.temp) %*% Cov.mat.temp + diag(n.subj * q))
	exp.bhat = full.inv.Cov.mat %*% t(Cov.mat.temp) %*% resid.vec.current
	
	gamma_physical <- matrix(0, nrow = n.subj, ncol = q)
	for(i in 1:n.subj) {
		b_i <- exp.bhat[(q*(i-1)+1):(q*i)]
		gamma_physical[i, ] <- temp.mat %*% b_i
	}
	
	fit = list()
	fit$fixed = beta.opt
	fit$stddev = sqrt(diag(Cov.Mat.RE))
	fit$BIC = BIC.values
	fit$t.frac = opt.frac
	fit$sigma.2 = sigma.2.opt
	fit$corr = round(diag(1/(fit$stddev+0.00000001))%*%Cov.Mat.RE%*%diag(1/(fit$stddev+0.00000001)),6)
	fit$bhat = gamma_physical
	
	return(fit)
}
