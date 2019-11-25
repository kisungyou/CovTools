#' One-Sample Covariance Test by Cai and Ma (2013)
#'
#' Given data, it performs 1-sample test for Covariance where
#' the null hypothesis is
#' \deqn{H_0 : \Sigma_n = \Sigma_0}
#' where \eqn{\Sigma_n} is the covariance of data model and \eqn{\Sigma_0} is a
#' hypothesized covariance based on a procedure proposed by Cai and Ma (2013).
#'
#' @param data an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param Sigma0 a \eqn{(p\times p)} given covariance matrix.
#' @param alpha level of significance.
#'
#' @return a named list containing: \describe{
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
#' CovTest1.2013Cai(data)
#' }
#'
#' @references
#' \insertRef{cai_optimal_2013}{CovTools}
#'
#' @export
CovTest1.2013Cai <- function(data, Sigma0=diag(ncol(data)), alpha=0.05){
  #---------------------------------------------------------------------------------
  ## PREPROCESSING ON INPUTS AND PARAMETERS
  # 1. valid data matrix
  if (!check_datamatrix(data)){
    stop("* CovTest1.2013Cai : an input matrix is invalid.")
  }
  n = nrow(data)
  p = ncol(data)
  if ((nrow(data)==1)||(ncol(data)==1)){stop("* CovTest1.2013Cai : invalid input matrix X.")}
  # 2. valid Sigma0
  if ((!check_sqmat(Sigma0))||(!check_pd(Sigma0))||(!isSymmetric(Sigma0,tol=sqrt(.Machine$double.eps)))||(nrow(Sigma0)!=p)){
    stop("* CovTest1.2013Cai : a given matrix for null hypothesis 'Sigma0' is invalid.")
  }
  # 3. alpha
  if ((length(alpha)!=1)||(alpha<=0)||(alpha>=1)){
    stop("* CovTest1.2013Cai : 'alpha' should be a real number in (0,1).")
  }

  #---------------------------------------------------------------------------------
  ## Adjust the matrix
  scaler = get_invroot(Sigma0)
  X.centered = scale(data, center=TRUE, scale=FALSE)
  X.adjusted = (matrix(X.centered,nrow=n) %*% scaler)

  #---------------------------------------------------------------------------------
  ## MAIN COMPUTATION
  output = test1.Cai13(X.adjusted, alpha)

  #---------------------------------------------------------------------------------
  ## RESULTS
  return(output)
}

# [2013.Cai] Optimal hypothesis testing for high-dimensional covar --------
# Tony Cai and Zongming Ma. (2013) Optimal hypothesis testing for high dimensional covariance matrices.
# If (Tn > Trej), reject Null Hypothesis (Sigma=I)
#' @keywords internal
#' @noRd
test1.Cai13 <- function(X, alpha){
  # 3. parameters and hXiXj
  n = nrow(X)
  p = ncol(X)
  hXiXj = rcpptest1_cai11(X)
  # 4. post-adjusting
  statistic = (hXiXj*2)/(n*(n-1)) # test statistic
  threshold = (qnorm(1-alpha))*2*sqrt((p*(p+1))/(n*(n-1)))
  reject = (statistic > threshold)
  # 5. return
  return(list(statistic=statistic, threshold=threshold, reject=reject))
}
