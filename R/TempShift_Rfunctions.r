KernelRegression <- function(x, y, sigma, kernelPar, k=GaussianKernel) {
	K <- k(x, x, kernelPar) # K <- k(x, x, 1,1)
	L <- chol(K + sigma^2*diag(x=1, nrow=nrow(K)))
	
	alpha <- backsolve(L, forwardsolve(t(L), y))
	predictNew <- function(newX) {
		if (is.vector(newX)) {
			newX <- t(as.matrix(newX))
		}
		# estimate mean
		newK <- k(x, newX, kernelPar)
		pred.mu <- t(newK) %*% alpha
		
		# estimate variance
		v = forwardsolve(t(L), newK)
		pred.var <- k(newX, newX, kernelPar)-t(v) %*% v
		
		return(list(Mean=pred.mu, Variance=pred.var))
	}

	logLik <- -1/2*t(y) %*% alpha - sum(log(diag(L)))-nrow(K)/2*log(2*pi)
	Fitted <- c(predictNew(x), LogLik=logLik)
	return(list(Fitted=Fitted, newFunction=predictNew))
}

selParametersByLikelihood <- function(x, y, initialPars, k=GaussianKernel) { # initialPars=c(l=<number>, sigmaF=<number>, sigma=<number>)
	calLogLik <- function(hp) {
		hyPars <- as.list(hp[names(hp) != "sigma"])
		assign("sigma", hp["sigma"])
		
		K <- k(x, x, kernelPar=hyPars)
		L <- chol(K + sigma^2*diag(x=1, nrow=nrow(K)))
		alpha <- backsolve(L, forwardsolve(t(L), y))
		logLik <- -1/2*t(y) %*% alpha - sum(log(diag(L)))-length(x)/2*log(2*pi)
		return(-logLik)
	}
	out <- optim(initialPars, calLogLik, method="BFGS")
	names(out$par) <- names(initialPars)
	return(out)
}

#####################################################################################
# multiple group models
multiGroupSelParametersShiftModelConstraint <- function(xLs, yLs, initialKernelPars, initialDTs, lowerDTs, upperDTs, k=GaussianKernel,
														lowerKernelPars=rep(-Inf, length(initialKernelPars)), upperKernelPars=rep(Inf, length(initialKernelPars))) {
	# initialPars=c(l=<number>, sigmaF=<number>, sigma=<number>) initialDTs=c(dT1=<number>, dT2=<number>, ...))
	names(initialDTs) <- paste("dT", 1:length(initialDTs), sep="")
	calLogLik <- function(hp) {
		hyPars <- as.list(hp[names(hp) != "sigma"])
		assign("sigma", hp["sigma"])
		dTLs <- hp[grep("dT", names(hp))]
		
		x <- do.call(c, lapply(1:length(xLs), function(i){xLs[[i]]-c(0, dTLs)[i]}))
		y <- do.call(c, yLs)
		K <- k(x, x, hyPars) # K <- k(x, x, 1,1)
		L <- chol(K + sigma^2*diag(x=1, nrow=nrow(K)))
	
		alpha <- backsolve(L, forwardsolve(t(L), y))
		logLik <- -1/2*t(y) %*% alpha - sum(log(diag(L)))-length(x)/2*log(2*pi)
		return(-logLik)
	}

	initialPars <- c(initialKernelPars, initialDTs)
	out <- optim(initialPars, calLogLik, method="L-BFGS-B", lower=c(lowerKernelPars, lowerDTs), upper=c(upperKernelPars, upperDTs))
	names(out$par) <- names(initialPars)
	return(out)
}

multiGroupSelCommonParameterIndependenceModel <- function(xLs, yLs, initialPars, k=GaussianKernel) { # initialPars=c(l=<number>, sigmaF=<number>, sigma=<number>)
	calLogLik <- function(hp) {
		hyPars <- as.list(hp[names(hp) != "sigma"])
		assign("sigma", hp["sigma"])
		
		x <- do.call(c, xLs)
		y <- do.call(c, yLs)
		KLs <- lapply(1:length(xLs), function(i){
			K <- k(xLs[[i]], xLs[[i]], hyPars)
			return(K)
		})
		
		groupN <- sapply(KLs, nrow)
		K <- matrix(0, nrow=length(x), ncol=length(x))
		for (i in 1:length(xLs)) {
			if (i==1){
				K[1:length(xLs[[1]]), 1:length(xLs[[1]])] <- KLs[[i]]
			} else {
				groupInd <- (sum(groupN[1:(i-1)])+1):sum(groupN[1:i])
				K[groupInd, groupInd] <- KLs[[i]]
			}	
		}
		
		L <- chol(K + sigma^2*diag(x=1, nrow=nrow(K)))
	
		alpha <- backsolve(L, forwardsolve(t(L), y))
		logLik <- -1/2*t(y) %*% alpha - sum(log(diag(L)))-length(x)/2*log(2*pi)
		return(-logLik)
	}

	out <- optim(initialPars, calLogLik, method="BFGS")
	names(out$par) <- names(initialPars)
	return(out)
}

fitShiftModel <- function(xLs, yLs, KernelPar0, lowerDTs, upperDTs, dTs0=rep(0, length(xLs)-1), k=GaussianKernel, lowerKernelPars=rep(-Inf, length(KernelPar0)), upperKernelPars=rep(Inf, length(KernelPar0))) { # initialPars must have names
	# shiftModel
	shiftPar <- multiGroupSelParametersShiftModelConstraint(xLs, yLs, initialKernelPars=KernelPar0, initialDTs=dTs0, lowerDTs, upperDTs, k,
															lowerKernelPars, upperKernelPars)
	dTLs <- shiftPar$par[grep("dT", names(shiftPar$par))]
	x <- do.call(c, lapply(1:length(xLs), function(i){xLs[[i]]-c(0, dTLs)[i]}))
	y <- do.call(c, yLs)
	shiftFit <- KernelRegression(x, y, sigma=shiftPar$par["sigma"], kernelPar=as.list(shiftPar$par), k)
	return(list(Par=shiftPar$par, Model=shiftFit))
}

fitIndependenceModel <- function(xLs, yLs, KernelPar0, k=GaussianKernel) { # initialPars must have names
	# independence model common par
	indPar <- multiGroupSelCommonParameterIndependenceModel(xLs, yLs, KernelPar0, k)
	fitLs <- lapply(1:length(xLs), function(i) {
		KernelRegression(xLs[[i]], yLs[[i]], sigma=indPar$par["sigma"], kernelPar=as.list(indPar$par), k)
	})
	names(fitLs) <- names(xLs)
	return(list(Par=indPar$par, Model=fitLs))
}

fitNoShiftModel_shiftPar <- function(xLs, yLs, ShiftPar, k=GaussianKernel) { # initialPars must have names
	# shiftModel
	ShiftPar <- as.list(ShiftPar)
	fit0 <- KernelRegression(unlist(xLs), unlist(yLs), sigma=ShiftPar$sigma, kernelPar=as.list(ShiftPar), k)

	return(list(Model=fit0))
}

fitAllModels <- function(xLs, yLs, KernelPar0, lowerDTs, upperDTs, dTs0=rep(0, length(xLs)-1), k=GaussianKernel, lowerKernelPars=rep(-Inf, length(KernelPar0)), upperKernelPars=rep(Inf, length(KernelPar0))) { # initialPars must have names
	# independence model common par
	indPar <- multiGroupSelCommonParameterIndependenceModel(xLs, yLs, KernelPar0, k)
	fitLs <- lapply(1:length(xLs), function(i) {
		KernelRegression(xLs[[i]], yLs[[i]], sigma=indPar$par["sigma"], kernelPar=as.list(indPar$par), k)
	})
	names(fitLs) <- names(xLs)
	
	# shiftModel
	shiftPar <- multiGroupSelParametersShiftModelConstraint(xLs, yLs, initialKernelPars=KernelPar0, initialDTs=dTs0, lowerDTs, upperDTs, k,
															lowerKernelPars, upperKernelPars)
	dTLs <- shiftPar$par[grep("dT", names(shiftPar$par))]
	x <- do.call(c, lapply(1:length(xLs), function(i){xLs[[i]]-c(0, dTLs)[i]}))
	y <- do.call(c, yLs)
	shiftFit <- KernelRegression(x, y, sigma=shiftPar$par["sigma"], kernelPar=as.list(shiftPar$par), k)
	
	# no-shift model
	fit0 <- KernelRegression(unlist(xLs), unlist(yLs), sigma=shiftPar$par["sigma"], kernelPar=as.list(shiftPar$par), k)
	
	return(list(IndependenceModel=list(Par=indPar$par, Model=fitLs),
				ShiftModel=list(Par=shiftPar$par, Model=shiftFit),
				NoShiftModel=list(Model=fit0)))
}

####################################
# GP simulator
library(mvtnorm)
GPgenerator <- function(n, x, kernelPars, sigma, k=GaussianKernel) {
	K <- k(x, x, kernelPars)
	t(mvtnorm::rmvnorm(n, mean=rep(0, length(x)), sigma=K+diag(sigma^2, length(x))))
}

