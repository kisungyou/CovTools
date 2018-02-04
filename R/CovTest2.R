#' Two-Sample Tests for Covariance Matrices
#'
#' Given two sets of data, \code{CovTest2} performs 2-sample test for Covariance where
#' the null hypothesis is
#' \deqn{H_0 : \Sigma_1 = \Sigma_2}
#' where \eqn{\Sigma_1} and \eqn{\Sigma_2} represent true (unknown) covariance
#' for each dataset.
#'
#' @param X an \code{(m-by-p)} matrix where each row is an observation from the first dataset.
#' @param Y an \code{(n-by-p)} matrix where each row is an observation from the second dataset.
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
#' ## generate 2 datasets from multivariate normal with identical covariance.
#' data1 = mvtnorm::rmvnorm(100, sigma=diag(5))
#' data2 = mvtnorm::rmvnorm(200, sigma=diag(5))
#'
#' ## run test
#' CovTest2(data1, data2)
#'
#'
#' @references [Cai13] Cai, T., Liu, W., and Xia, Y. (2013) \emph{Two-Sample Covariance Matrix Testing and Support Recovery in
#' High-Dimensional and Sparse Settings.} Journal of American Statistical Association, Vol.108(501):265-277.
#' @export
CovTest2 <- function(X, Y, alpha=0.05, method=c("Cai13")){
  ## PREPROCESSING ON INPUTS AND PARAMETERS
  # 1. valid data matrix
  if (!check_datamatrix(X)){stop("* CovTest2 : an input data matrix X is invalid.")}
  if (!check_datamatrix(Y)){stop("* CovTest2 : an input data matrix Y is invalid.")}
  if ((nrow(X)==1)||(ncol(X)==1)){stop("* CovTest2 : invalid input matrix X.")}
  if ((nrow(Y)==1)||(ncol(Y)==1)){stop("* CovTest2 : invalid input matrix Y.")}
  if (ncol(X)!=ncol(Y)){stop("* CovTest2 : two inpu matrices X and Y should have same number of columns.")}
  # 2. alpha
  if ((length(alpha)!=1)||(alpha<=0)||(alpha>=1)){stop("* CovTest2 : 'alpha' should be a real number in (0,1).")}
  # 2. method : THIS SHOULD BE UPDATED EVERYTIME A METHOD IS ADDED
  if (missing(method)){method="Cai13"}
  method = match.arg(method)

  ## MAIN COMPUTATION
  output = switch(method,
                  Cai13 = test2.Cai13(X, Y, alpha)
  )

  ## RETURN OUTPUT
  return(output)
}
