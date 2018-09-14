# [2018.Lee] Maximum Pairwise Bayes Factors for Covariance Structure Testing
#' @keywords internal
#' @noRd
bayestest1.Lee18 <- function(X, params){
  #################################################################
  ## Params
  parnames = names(params)
  if ("a0" %in% parnames){    a0 = params$a0  } else {    a0 = 2.0} ## a0
  if ("b0" %in% parnames){    b0 = params$b0  } else {    b0 = 2.0} ## a0
  if ("gamma" %in% parnames){gamma= params$gamma}else{ gamma = 1.0}
  ## Parameter Value Warning
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
      log.BF.mat[i,j] = a0*log(b0) - lgamma(a0) +
        1/2 * log(gamma/(1+gamma)) + lgamma(n/2 + a0) +
        1/2 * sum((Xi)^2) - (n/2 + a0) * log(1/2 * ( sum((Xi)^2) - sum(Xi*Xj)^2/sum((Xj)^2) /(1+gamma) ) + b0)
    }
  }
  diag(log.BF.mat) = -Inf # just to fill out the diagonal parts
  output = list()
  output$log.BF.mat = log.BF.mat
  return(output)
}
