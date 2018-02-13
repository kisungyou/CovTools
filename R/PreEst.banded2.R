# Original name : Lee17
#' Bayesian Estimation of Precision Matrix with a Banding Structure (Lee 2017)
#'
#' @param X an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param upperK upper bound of bandwidth \eqn{k}
#' @param logpi log of prior distribution for bandwidth \eqn{k}. Default is a function proportional to \eqn{-k^4}.
#'
#' @return a named list containing: \describe{
#' \item{C}{a \eqn{(p\times p)} estimated precision matrix.}
#' }
#'
#' @examples
#' ## generate data from multivariate normal with Identity precision.
#' data = mvtnorm::rmvnorm(100, sigma=diag(10))
#'
#' ## compare different K
#' out1 <- PreEst.banded2(data, upperK=1)
#' out2 <- PreEst.banded2(data, upperK=3)
#' out3 <- PreEst.banded2(data, upperK=5)
#'
#' ## visualize
#' par(mfrow=c(2,2), pty="s")
#' image(pracma::flipud(diag(10)),main="Original Precision")
#' image(pracma::flipud(out1$C), main="banded2::upperK=1")
#' image(pracma::flipud(out2$C), main="banded2::upperK=3")
#' image(pracma::flipud(out3$C), main="banded2::upperK=5")
#'
#' @references
#' \insertRef{lee_estimating_2017}{CovTools}
#'
#' @rdname PreEst.banded2
#' @export
PreEst.banded2 <- function(X, upperK=floor(ncol(X)/2), logpi=function(k){-k^4}){
  #-----------------------------------------------------
  ## PREPROCESSING
  #   1. typecheck : data
  fname    = "PreEst.banded2"
  checker1 = invisible_datamatrix(X, fname)
  #   2. upperK : upper bound of bandwidth
  Lee17.upperK = as.integer(upperK)
  if ((length(Lee17.upperK)!=1)||(Lee17.upperK<1)||(Lee17.upperK>ncol(X))||(abs(Lee17.upperK-round(Lee17.upperK))>sqrt(.Machine$double.eps))){
    stop("* PreEst.banded2 : 'upperK' should be an integer in [1,ncol(X)].")
  }
  #   3. logpi : log of prior distribution for bandwidth k.
  Lee17.logpi = logpi
  if (!is.function(Lee17.logpi)){
    stop("* PreEst.banded2 : 'logpi' should be a function.")
  }

  #-----------------------------------------------------
  ## PREPROCESSING
  output      = preest.Lee17(Lee17.upperK, X, logpi=Lee17.logpi)

  #-----------------------------------------------------
  ## RETURN
  return(output)
}
