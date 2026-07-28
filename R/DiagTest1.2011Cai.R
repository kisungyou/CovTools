#' One-Sample Diagonality Test by Cai and Jiang (2011)
#'
#' Given data, it performs 1-sample test for diagonal entries of a Covariance matrix where
#' the null hypothesis is
#' \deqn{H_0 : \sigma_{ij} = 0~\mathrm{for any}~i \neq j}
#' and alternative hypothesis is \eqn{H_1 : ~\mathrm{not}~H_0}
#' with \eqn{\Sigma_n = (\sigma_{ij})} based on a procedure proposed by Cai and Jiang (2011).
#' Each column is centered before the sample correlations are computed, so the
#' test is invariant to translations of the observations.
#'
#' @param data an \eqn{(n\times p)} data matrix where each row is an observation.
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
#' data = matrix(rnorm(100*pdim), ncol=pdim)
#'
#' ## run test with different alpha values
#' DiagTest1.2011Cai(data, alpha=0.01)
#' DiagTest1.2011Cai(data, alpha=0.05)
#' DiagTest1.2011Cai(data, alpha=0.10)
#' }
#'
#' @references
#' \insertRef{cai_limiting_2011}{CovTools}
#'
#' @export
DiagTest1.2011Cai <- function(data, alpha=0.05){
  ###########################################################################
  # Preprocessing : Inputs
  # 1. data
  if (!check_datamatrix(data)){
    stop("* DiagTest1.2011Cai : an input matrix is invalid.")
  }
  n = nrow(data)
  p = ncol(data)
  if ((nrow(data)==1)||(ncol(data)==1)){stop("* DiagTest1.2011Cai : invalid input matrix X.")}
  # 2. alpha
  if ((length(alpha)!=1)||(alpha<=0)||(alpha>=1)){
    stop("* DiagTest1.2011Cai : 'alpha' should be a real number in (0,1).")
  }

  ###########################################################################
  ## MAIN COMPUTATION
  output = diagtest1.Cai11(data, alpha)

  ###########################################################################
  ## RETURN OUTPUT
  return(output)
}


# auxiliary function ------------------------------------------------------
## Cai and Jiang (2011)'s test
## LIMITING LAWS OF COHERENCE OF RANDOM MATRICES WITH APPLICATIONS TO TESTING COVARIANCE STRUCTURE AND CONSTRUCTION OF COMPRESSED SENSING MATRICES
#' @keywords internal
#' @noRd
diagtest1.Cai11 <- function(X, alpha){
  p = ncol(X)
  n = nrow(X)

  X = sweep(X, 2L, colMeans(X), FUN="-", check.margin=FALSE)
  column.ss = colSums(X^2)
  if ((any(!is.finite(column.ss)))||(any(column.ss<=0))){
    stop("* DiagTest1.2011Cai : every column should have positive sample variation.")
  }
  Rho = crossprod(X)/sqrt(outer(column.ss, column.ss))
  diag(Rho) = 0
  Ln = max(abs(Rho))
  Tn = n*(Ln)^2 - 4*log(p) + log(log(p))

  res = list()
  res$statistic = Tn
  res$threshold = -log(8*pi) - 2*log(log(1/(1-alpha)))
  res$reject = (Tn > res$threshold)
  return(res)
}
