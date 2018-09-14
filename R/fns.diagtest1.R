## Lan et al. (2015)'s test for diagonality
## Testing the Diagonality of a Large Covariance Matrix in a Regression Setting
#' @keywords internal
#' @noRd
diagtest1.Lan15 <- function(X, alpha){
  p = ncol(X)
  n = nrow(X)
  z.val = qnorm(1 - alpha)

  Sg.hat = 0
  for(i in 1:n){
    X.i = matrix(X[i,], ncol=1)
    Sg.hat = Sg.hat + X.i%*%t(X.i)/n
  }
  Sg.hat2 = Sg.hat^2

  Bias.hat = sqrt(n)/(2*p^{3/2}) * ( sum(diag(Sg.hat))^2 - sum(diag(Sg.hat2)) )
  M2p.hat = n/(p*(n+2)) * sum(diag(Sg.hat2))
  Tstat = 0
  for(j1 in 1:(p-1)){
    for(j2 in (j1+1):p){
      Tstat = Tstat + (n/p)^{3/2} * Sg.hat2[j1, j2]
    }
  }
  Tn = (Tstat - Bias.hat) / (sqrt(n/p)*M2p.hat)

  res = list()
  res$statistic = Tn
  res$threshold = z.val
  res$reject = ( Tn > z.val )
  return( res )
}

## Cai and Jiang (2011)'s test
## LIMITING LAWS OF COHERENCE OF RANDOM MATRICES WITH APPLICATIONS TO TESTING COVARIANCE STRUCTURE AND CONSTRUCTION OF COMPRESSED SENSING MATRICES
#' @keywords internal
#' @noRd
diagtest1.Cai11 <- function(X, alpha){
  p = ncol(X)
  n = nrow(X)

  Rho = matrix(0, p,p)
  for(i in 1:p){
    X.i = matrix(X[,i], ncol=1)
    for(j in (1:p)[-i]){
      X.j = matrix(X[,j], ncol=1)
      Rho[i,j] = sum(X.i*X.j) / sqrt( sum( X.i^2 )*sum( X.j^2 ) )
    }
  }
  Ln = max(abs(Rho))
  Tn = n*(Ln)^2 - 4*log(p) + log(log(p))

  res = list()
  res$statistic = Tn
  res$threshold = -log(8*pi) - 2*log(log(1/(1-alpha)))
  res$test = (Tn > res$threshold)
  return(res)
}

# - Kyoungjae Lee, Lizhen Lin and David Dunson. (2018) Maximum Pairwise Bayes Factors for Covariance Structure Testing. [https://arxiv.org/abs/1809.03105]
#' @keywords internal
#' @noRd
diagtest1.mxPBF <- function(X, params){
  #################################################################
  # Params
  parnames = names(params)
  if ("a0" %in% parnames){    a0 = params$a0  } else {    a0 = 2.0} ## a0
  if ("b0" %in% parnames){    b0 = params$b0  } else {    b0 = 2.0} ## a0
  if ("gamma" %in% parnames){gamma= params$gamma}else{ gamma = 1.0}

  if ((length(a0)!=1)||(a0<=0)){stop("* mxPBF : 'a0' should be nonnegative number.")}
  if ((length(b0)!=1)||(b0<=0)){stop("* mxPBF : 'b0' should be nonnegative number.")}
  if ((length(gamma)!=1)||(gamma<=0)){stop("* mxPBF : 'gamma' should be nonnegative number.")}

  #################################################################
  ## MAIN RUN BY JAY
  p = ncol(X)
  n = nrow(X)
  log.BF.mat = matrix(0, ncol=p,nrow=p) # log Bayes factors

  for(i in 1:p){
    Xi = matrix(X[,i], ncol=1)
    for(j in (1:p)[-i]){
      Xj = matrix(X[,j], ncol=1)
      log.BF.mat[i,j] = 1/2 * log(gamma/(1+gamma)) -
        (n/2 + a0) * ( log(( sum((Xi)^2) - sum(Xi*Xj)^2/sum((Xj)^2) /(1+gamma) ) + 2*b0) - log( sum((Xi)^2) + 2*b0) )
    }
  }
  diag(log.BF.mat) = -Inf # just to fill out the diagonal parts
  return(log.BF.mat)
}
