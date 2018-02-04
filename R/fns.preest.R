# [Yuan07] Penalized Likelihood type of L1 --------------------------------
#' @keywords internal
#' @noRd
preest.Yuan07.once.matrix <- function(X, lambda){
  S    = cov(X)
  Chat = rcpp_ADMMprecision(S, lambda);

  output = list()
  output$C = matrix(Chat, nrow=nrow(S))
  return(output)
}
#' @keywords internal
#' @noRd
preest.Yuan07.once.BIC <- function(X, lambda){
  n    = nrow(X)
  p    = ncol(X)
  S    = cov(X)
  C    = rcpp_ADMMprecision(S, lambda);

  term1 = -log(det(C))
  term2 = sum(diag(C%*%S))
  term3 = 0
  sqrtprecision = sqrt(.Machine$double.eps)
  for (i in 1:p){
    for (j in i:p){
      if (abs(C[i,j]) <= sqrtprecision){
        term3 = term3+1.0
      }
    }
  }
  BICscore = -term1+term2+((log(n)/n)*term3)
  return(BICscore)
}
#' @keywords internal
#' @noRd
preest.Yuan07.plgrid <-function(X, plgrid, nCore){
  S = cov(X)
  if (nCore==1){
    cl = makeCluster(1)
    registerDoParallel(cl)
  } else {
    cl = makeCluster(nCore)
    registerDoParallel(cl)
  }
  # parallel computation let's go
  # I had a similar problem and I solved it by adding .noexport = c(<Functions that were implemented in C++>) to the foreach.
  itforeach=NULL
  BICs = foreach (itforeach=1:length(plgrid), .combine=rbind) %dopar% {
    preest.Yuan07.once.BIC(X,plgrid[itforeach])
  }
  # stop cluster
  stopCluster(cl)
  # find which should be the best one
  idxmin = which.min(BICs)
  if (length(idxmin)>1){idxmin = idxmin[1]}
  # optimal threshold value
  lambdaopt = plgrid[idxmin]
  Copt = preest.Yuan07.once.matrix(X, lambdaopt)

  # return output
  output = list()
  output$C = matrix(Copt$C, nrow=nrow(S))
  output$BIC = data.frame(lambda=plgrid, BIC=BICs)
  rownames(output$BIC) = NULL
  return(output)
}

# [Banerjee06] use heuristic for lambda determination ---------------------
#' @keywords internal
#' @noRd
preest.Banerjee06 <- function(X,confidence){
  # we need \rho value
  n = nrow(X)
  S = cov(X); diagS = (diag(S)); SS = outer(diagS,diagS);
  maxSS = max(SS);
  tval  = qt(0.5+(confidence/2),df=(n-2))
  rho   = tval*maxSS/sqrt(n-2+(tval^2))

  # run ADMM
  output = preest.Yuan07.once.matrix(X, rho)
  return(output)
}


# [Lee.17] approximated posterior mean of banded precision matrix  --------
#######################################################################
# approximated posterior mean of banded precision matrix based on Chol. decomp. (Lee and Lee, 2017+)
# K : upper bound of bandwidth
# X : n-by-p data matrix
# logpi : log of prior distribution for bandwidth k. Default is a function proportional to -k^4.
#
# refer : Estimating Large Precision Matrices via Modified Cholesky Decomposition (https://arxiv.org/abs/1707.01143#)
#######################################################################
# main function
#' @keywords internal
#' @noRd
preest.Lee17<- function(K, X, logpi=function(k){ -k^4 }){
  X = scale(X, center = TRUE, scale =FALSE)

  n = nrow(X); p = ncol(X)
  A = matrix(0, p, p)
  Dinv = diag(1, nrow=p)

  logprob = rep(0, K)
  for(k in 1:K){
    logprob[k] = Lee17.logpostpi(k, X, logpi)
  }
  khatLL = which(logprob == max(logprob)) # posterior mode of k

  dhat1 = Lee17.dhat(1, khatLL, X)
  Dinv[1, 1] = ((n-2)/(n*dhat1))
  for(j in 2:p){
    nj = n - min(j-1, khatLL) - 2
    ind = c(max(1,j-khatLL):(j-1))

    Xj = X[, j]
    Zj = Lee17.Z(j, khatLL, X)
    Covhat = t(Zj)%*%matrix(Xj)/n
    dhatj = Lee17.dhat(j, khatLL, X)

    A[j, ind] = solve(Lee17.VhatZ(j, khatLL, X))%*%Covhat
    Dinv[j, j] = ((nj)/(n*dhatj))
  }
  OmegaLL = (diag(1, nrow=p) - t(A))%*%Dinv%*%(diag(1, nrow=p) - A)

  output = list()
  output$C = OmegaLL
  return(output)
}
#######################################################################
# auxiliary functions
#######################################################################
# Z_j(k) : n X k matrix
#' @keywords internal
#' @noRd
Lee17.Z <- function(j, k, X){
  ind = c(max(1, j-k):(j-1))
  res = X[, ind]
  if(j == 2) res = matrix(res)
  return(res)
}
# \hat{Var}(Z_j(k)) : k X k matrix
#' @keywords internal
#' @noRd
Lee17.VhatZ <- function(j, k, X){
  Zj = Lee17.Z(j, k, X)
  res = t(Zj)%*%Zj/nrow(X)
  return(res)
}
# \hat{d}_{jk}
#' @keywords internal
#' @noRd
Lee17.dhat <- function(j, k, X){
  Xj = X[, j]
  if(j == 1){
    res = sum(Xj^2)/nrow(X)   # scalar
  }else{
    Zj = Lee17.Z(j, k, X)
    VhatX = sum(Xj^2)/nrow(X)   # scalar
    Covhat = t(Zj)%*%matrix(Xj)/nrow(X)   # k X 1 matrix
    res = VhatX - t(Covhat)%*%solve(Lee17.VhatZ(j, k, X))%*%Covhat
  }
  return(res)
}
# log of posterior probability for k : \log(\pi(k | \bfX_n)) (Lee and Lee, 2016)
#' @keywords internal
#' @noRd
Lee17.logpostpi <- function(k, X, logpi){
  p = ncol(X); n = nrow(X)
  res = logpi(k)
  for(j in 2:p){
    nj = n - min(j-1,k) - 2
    res = res - 0.5*determinant(n*Lee17.VhatZ(j, k, X)/(2*pi), logarithm=T)$modulus + lgamma(0.5*nj) - 0.5*nj*log(0.5*n*Lee17.dhat(j, k, X))
  }
  return(res)
}


# [Banerjee14] Bayesian estimation for banded precision matrix bas --------
#######################################################################
# Bayesian estimation for banded precision matrix based on G-Wishart prior (Banerjee and Ghosal, 2014)
# K : upper bound of bandwidth
# X : n-by-p data matrix
# delta : hyperparameter for G-Wishart prior. Default value is 10. It has to be larger than 2.
# logpi : log of prior distribution for bandwidth k. Default is a function proportional to -k^4.
#
# refer : Banerjee and Ghosal, 2014. Posterior convergence rates for estimating large precision matrices using graphical models
#######################################################################
# branching case for Jay's code
#' @keywords internal
#' @noRd
preest.Banerjee14 <- function(K, X, delta, logpi, type){
  if (type=="Stein"){
    outC = preest.Banerjee14.GWishartE1(K, X, delta, logpi)
  } else if (type=="Squared"){
    outC = preest.Banerjee14.GWishartE2(K, X, delta, logpi)
  }

  output = list()
  output$C = outC
  return(output)
}
# main functions
# (1) Bayes estimator under the Stein's loss
#' @keywords internal
#' @noRd
preest.Banerjee14.GWishartE1 <- function(K, X, delta, logpi){
  X = scale(X, center = TRUE, scale =FALSE)

  logprobBG = rep(0, K)
  for(k in 1:K){
    logprobBG[k] = preest.Banerjee14.logJG(k, X, delta) + logpi(k)
  }
  khatBG = which(logprobBG == max(logprobBG)) # choice of bandwidth

  n = nrow(X); p = ncol(X)
  S = t(X)%*%X/n

  OmegaBG1 = matrix(0, p,p)
  OmegaBG1[1:(khatBG+1), 1:(khatBG+1)] = solve(diag(1/n, khatBG+1) + S[1:(khatBG+1), 1:(khatBG+1)])
  for(j in 2:(p-k)){
    OmegaBG1[j:(j+khatBG), j:(j+khatBG)] = OmegaBG1[j:(j+khatBG), j:(j+khatBG)] + solve(diag(1/n, khatBG+1) + S[j:(j+ khatBG), j:(j+khatBG)])
    OmegaBG1[j:(j+khatBG-1), j:(j+khatBG-1)] = OmegaBG1[j:(j+khatBG-1), j:(j+khatBG-1)] - solve(diag(1/n,khatBG) + S[j:(j+khatBG-1), j:(j+khatBG-1)])
  }
  OmegaBG1 = OmegaBG1*(delta+n-2)/n

  return(OmegaBG1)
}
# (2) Bayes estimator uner the Squared-error loss
#' @keywords internal
#' @noRd
preest.Banerjee14.GWishartE2 <- function(K, X, delta, logpi){
  X = scale(X, center = TRUE, scale =FALSE)

  logprobBG = rep(0, K)
  for(k in 1:K){
    logprobBG[k] = preest.Banerjee14.logJG(k, X, delta) + logpi(k)
  }
  khatBG = which(logprobBG == max(logprobBG)) # choice of bandwidth

  n = nrow(X); p = ncol(X)
  S = t(X)%*%X/n

  OmegaBG2 = matrix(0, p,p)
  OmegaBG2[1:(khatBG+1), 1:(khatBG+1)] = solve(diag(1/n,khatBG+1) + S[1:(khatBG+1), 1:(khatBG+1)])*(delta+khatBG+n)/n
  for(j in 2:(p-khatBG)){
    OmegaBG2[j:(j+khatBG), j:(j+khatBG)] = OmegaBG2[j:(j+khatBG), j:(j+khatBG)] + solve(diag(1/n,khatBG+1) + S[j:(j+khatBG), j:(j+khatBG)])*(delta+khatBG+n)/n
    OmegaBG2[j:(j+khatBG-1), j:(j+khatBG-1)] = OmegaBG2[j:(j+khatBG-1), j:(j+khatBG-1)] - solve(diag(1/n,khatBG) + S[j:(j+khatBG-1), j:(j+khatBG-1)])*(delta+khatBG+n-1)/n
  }

  return(OmegaBG2)
}
########################################################################
# auxiliary functions
########################################################################
# logarithm of J_{G^k}(delta, n, I_p, nS) in Banerjee and Ghosal (2014) (p2123: eq (5.7))
#' @keywords internal
#' @noRd
preest.Banerjee14.logJG <- function(k, X, delta){
  n = nrow(X); p = ncol(X)
  S = t(X)%*%X/n

  res = 0
  res = res + sum(lgamma(0.5*(delta+n+(0:k))) - lgamma(0.5*(delta+(0:k))))
  res = res + (p-k-1)*(lgamma(0.5*(delta+n+k)) - lgamma(0.5*(delta+k)))
  SC1 = S[1:(1+k), 1:(1+k)]
  res = res - 0.5*(delta+n+k)*determinant(diag(1,nrow=k+1)+n*SC1, logarithm = T)$modulus
  for(j in 2:(p-k)){
    SSj = S[j:(j+k-1), j:(j+k-1)]
    SCj = S[j:(j+k), j:(j+k)]
    res = res + 0.5*(delta+n+k-1)*determinant(diag(1,nrow=k)+n*SSj, logarithm = T)$modulus - 0.5*(delta+n+k)*determinant(diag(1,nrow=k+1)+n*SCj, logarithm = T)$modulus
  }

  return(res)
}
