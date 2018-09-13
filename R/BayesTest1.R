#' Bayesian One-Sample Tests for Covariance Matrix
#'
#' @export
BayesTest1 <- function(data, Sigma0=diag(ncol(data)), method=c("mxPBF"), ...){
  ###########################################################################
  # Preprocessing : Inputs
  # 1. data
  if (!check_datamatrix(data)){
    stop("* BayesTest1 : an input matrix is invalid.")
  }
  n = nrow(data)
  p = ncol(data)
  if ((nrow(data)==1)||(ncol(data)==1)){stop("* BayesTest1 : invalid input matrix X.")}
  # 2. valid Sigma0
  if ((!check_sqmat(Sigma0))||(!check_pd(Sigma0))||(!isSymmetric(Sigma0,tol=sqrt(.Machine$double.eps)))||(nrow(Sigma0)!=p)){
    stop("* BayesTest1 : a given matrix for null hypothess 'Sigma0' is invalid.")
  }
  # 3. method
  method = match.arg(method)
  # 4. extra arguments
  extra.args = (list(...))

  ###########################################################################
  # Preprocessing : Adjust the data for testing I
  scaler = get_invroot(Sigma0)
  X.centered = scale(data, center=TRUE, scale=FALSE)
  X.adjusted = (matrix(X.centered,nrow=n) %*% scaler)

  ###########################################################################
  # Main Computation
  output = switch(method,
                  mxPBF = bayestest1.Lee18(X.adjusted, extra.args))

  ###########################################################################
  return(output)
}
