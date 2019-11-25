#' One-Sample Diagonality Test by Lan et al. (2015)
#'
#' Given data, it performs 1-sample test for diagonal entries of a Covariance matrix where
#' the null hypothesis is
#' \deqn{H_0 : \sigma_{ij} = 0~\mathrm{for any}~i \neq j}
#' and alternative hypothesis is \eqn{H_1 : ~\mathrm{not}~H_0}
#' with \eqn{\Sigma_n = (\sigma_{ij})} based on a procedure proposed by Lan et al. (2015).
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
#' DiagTest1.2015Lan(data, alpha=0.01)
#' DiagTest1.2015Lan(data, alpha=0.05)
#' DiagTest1.2015Lan(data, alpha=0.10)
#' }
#'
#' @references
#' \insertRef{lan_testing_2015}{CovTools}
#'
#' @export
DiagTest1.2015Lan <- function(data, alpha=0.05){
  ###########################################################################
  # Preprocessing : Inputs
  # 1. data
  if (!check_datamatrix(data)){
    stop("* DiagTest1.2015Lan : an input matrix is invalid.")
  }
  n = nrow(data)
  p = ncol(data)
  if ((nrow(data)==1)||(ncol(data)==1)){stop("* DiagTest1.2015Lan : invalid input matrix X.")}
  # 2. alpha
  if ((length(alpha)!=1)||(alpha<=0)||(alpha>=1)){
    stop("* DiagTest1.2015Lan : 'alpha' should be a real number in (0,1).")
  }

  ###########################################################################
  ## MAIN COMPUTATION
  output = diagtest1.Lan15(data, alpha)

  ###########################################################################
  ## RETURN OUTPUT
  return(output)
}



# auxiliary function ------------------------------------------------------
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
