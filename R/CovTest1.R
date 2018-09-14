#' One-Sample Tests for Covariance Matrices
#'
#' Given data, \code{CovTest1} performs 1-sample test for Covariance where
#' the null hypothesis is
#' \deqn{H_0 : \Sigma_n = \Sigma_0}
#' where \eqn{\Sigma_n} is the covariance of data model and \eqn{\Sigma_0} is a
#' hypothesized covariance.
#'
#' @param data an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param Sigma0 a \eqn{(p\times p)} given covariance matrix.
#' @param alpha level of significance.
#' @param method a name of test.
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
#' data = matrix(rnorm(100*5), nrow=100)
#'
#' ## run test
#' CovTest1(data, method="Cai13")
#' CovTest1(data, method="Srivastava14")
#' }
#'
#' @references
#' \insertRef{cai_optimal_2013}{CovTools}
#'
#' \insertRef{srivastava_tests_2014}{CovTools}
#'
#' @export
CovTest1 <- function(data, Sigma0=diag(ncol(data)), alpha=0.05, method=c("Cai13","Srivastava14")){
  #---------------------------------------------------------------------------------
  ## PREPROCESSING ON INPUTS AND PARAMETERS
  # 1. valid data matrix
  if (!check_datamatrix(data)){
    stop("* CovTest1 : an input matrix is invalid.")
  }
  n = nrow(data)
  p = ncol(data)
  if ((nrow(data)==1)||(ncol(data)==1)){stop("* CovTest1 : invalid input matrix X.")}
  # 2. valid Sigma0
  if ((!check_sqmat(Sigma0))||(!check_pd(Sigma0))||(!isSymmetric(Sigma0,tol=sqrt(.Machine$double.eps)))||(nrow(Sigma0)!=p)){
    stop("* CovTest1 : a given matrix for null hypothesis 'Sigma0' is invalid.")
  }
  # 3. alpha
  if ((length(alpha)!=1)||(alpha<=0)||(alpha>=1)){
    stop("* CovTest1 : 'alpha' should be a real number in (0,1).")
  }
  # 4. method : THIS SHOULD BE UPDATED EVERYTIME A METHOD IS ADDED
  if (missing(method)){method = "Cai13"}
  method = match.arg(method)

  #---------------------------------------------------------------------------------
  ## Adjust the matrix
  scaler = get_invroot(Sigma0)
  X.centered = scale(data, center=TRUE, scale=FALSE)
  X.adjusted = (matrix(X.centered,nrow=n) %*% scaler)

  #---------------------------------------------------------------------------------
  ## MAIN COMPUTATION
  output = switch(method,
                  Cai13 = test1.Cai13(X.adjusted, alpha),
                  Srivastava14 = test1.Srivastava14(X.adjusted, alpha)
  )

  #---------------------------------------------------------------------------------
  ## RESULTS
  return(output)
}
