#' One-Sample Covariance Test by Srivastava, Yanagihara, and Kubokawa (2014)
#'
#' Given data, it performs 1-sample test for Covariance where
#' the null hypothesis is
#' \deqn{H_0 : \Sigma_n = \Sigma_0}
#' where \eqn{\Sigma_n} is the covariance of data model and \eqn{\Sigma_0} is a
#' hypothesized covariance based on a procedure proposed by Srivastava, Yanagihara, and Kubokawa (2014).
#'
#' @param data an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param Sigma0 a \eqn{(p\times p)} given covariance matrix.
#' @param alpha level of significance.
#'
#' @return a named list containing \describe{
#' \item{statistic}{a test statistic value.}
#' \item{threshold}{rejection criterion to be compared against test statistic.}
#' \item{reject}{a logical; \code{TRUE} to reject null hypothesis, \code{FALSE} otherwise.}
#' }
#'
#' @examples
#' \dontrun{
#' ## generate data from multivariate normal with trivial covariance.
#' pdim = 5
#' data = matrix(rnorm(10*pdim), ncol=pdim)
#'
#' ## run the test
#' CovTest1.2014Srivastava(data)
#' }
#'
#' @references
#' \insertRef{srivastava_tests_2014}{CovTools}
#'
#' @export
CovTest1.2014Srivastava <- function(data, Sigma0=diag(ncol(data)), alpha=0.05){
  #---------------------------------------------------------------------------------
  ## PREPROCESSING ON INPUTS AND PARAMETERS
  # 1. valid data matrix
  if (!check_datamatrix(data)){
    stop("* CovTest1.2014Srivastava : an input matrix is invalid.")
  }
  n = nrow(data)
  p = ncol(data)
  if ((nrow(data)==1)||(ncol(data)==1)){stop("* CovTest1.2014Srivastava : invalid input matrix X.")}
  # 2. valid Sigma0
  if ((!check_sqmat(Sigma0))||(!check_pd(Sigma0))||(!isSymmetric(Sigma0,tol=sqrt(.Machine$double.eps)))||(nrow(Sigma0)!=p)){
    stop("* CovTest1.2014Srivastava : a given matrix for null hypothesis 'Sigma0' is invalid.")
  }
  # 3. alpha
  if ((length(alpha)!=1)||(alpha<=0)||(alpha>=1)){
    stop("* CovTest1.2014Srivastava : 'alpha' should be a real number in (0,1).")
  }

  #---------------------------------------------------------------------------------
  ## Adjust the matrix
  scaler = get_invroot(Sigma0)
  X.centered = scale(data, center=TRUE, scale=FALSE)
  X.adjusted = (matrix(X.centered,nrow=n) %*% scaler)

  #---------------------------------------------------------------------------------
  ## MAIN COMPUTATION
  output = test1.Srivastava14(X.adjusted, alpha)

  #---------------------------------------------------------------------------------
  ## RESULTS
  return(output)
}

# [2014.Srivastava] Tests for covariance matrices in high dimension with less sample size.
#' @keywords internal
#' @noRd
test1.Srivastava14 <- function(X, alpha){
  Sigma0 = diag(ncol(X))
  p = ncol(X)
  n = nrow(X)
  z.val = qnorm(1 - alpha)
  bar.X = matrix(colMeans(X), nrow = p, ncol = 1)
  Sn = cov(X)
  Yn = X - matrix(c(bar.X), nrow=n, ncol=p, byrow=T)
  a1.hat = sum(diag(Sn)) / p
  a2.hat = ( (n-1)*(n-2)*sum(diag( (n-1)^2*Sn%*%Sn )) - n*(n-1)*sum((rowSums(Yn^2))^2) + (sum(diag( (n-1)*Sn )))^2 ) /( p*n*(n-1)*(n-2)*(n-3) )
  T2 = (n-1)/2 * ( a2.hat - 2*a1.hat + 1 )

  res = list()
  res$statistic = T2
  res$threshold = z.val
  res$reject = ( T2 > z.val )
  return( res )
}
