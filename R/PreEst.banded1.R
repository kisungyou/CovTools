# Originla name : Banerjee14
#' Bayesian Estimation of a Banded Precision Matrix (Banerjee 2014)
#'
#' \code{PreEst.banded1} returns a Bayes estimator of the banded precision matrix using G-Wishart prior.
#' Stein’s loss or squared error loss function is used depending on the “loss” argument in the function.
#' The bandwidth is set at the mode of marginal posterior for the bandwidth parameter.
#'
#' @param X an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param upperK upper bound of bandwidth \eqn{k}.
#' @param delta hyperparameter for G-Wishart prior. Default value is 10. It has to be larger than 2.
#' @param logpi log of prior distribution for bandwidth \eqn{k}. Default is a function proportional to \eqn{-k^4}.
#' @param loss type of loss; either \code{"Stein"} or \code{"Squared"}.
#'
#' @return a named list containing: \describe{
#' \item{C}{a \eqn{(p\times p)} MAP estimate for precision matrix.}
#' }
#'
#' @examples
#' ## generate data from multivariate normal with Identity precision.
#' data = mvtnorm::rmvnorm(100, sigma=diag(10))
#'
#' ## compare different K
#' out1 <- PreEst.banded1(data, upperK=1)
#' out2 <- PreEst.banded1(data, upperK=3)
#' out3 <- PreEst.banded1(data, upperK=5)
#'
#' ## visualize
#' par(mfrow=c(2,2), pty="s")
#' image(pracma::flipud(diag(10)),main="Original Precision")
#' image(pracma::flipud(out1$C), main="banded1::upperK=1")
#' image(pracma::flipud(out2$C), main="banded1::upperK=3")
#' image(pracma::flipud(out3$C), main="banded1::upperK=5")
#'
#' @references
#' \insertRef{banerjee_posterior_2014}{CovTools}
#'
#' @rdname PreEst.banded1
#' @export
PreEst.banded1 <- function(X, upperK=floor(ncol(X)/2), delta=10.0, logpi=function(k){-k^4}, loss=c("Stein","Squared")){
  #-----------------------------------------------------
  ## PREPROCESSING
  #   1. typecheck : data
  fname    = "PreEst.banded1"
  checker1 = invisible_datamatrix(X, fname)
  #   2. upperK : upper bound of bandwidth
  Banerjee14.upperK = as.integer(upperK)
  if ((length(Banerjee14.upperK)!=1)||(Banerjee14.upperK<1)||(Banerjee14.upperK>ncol(X))||(abs(Banerjee14.upperK-round(Banerjee14.upperK))>sqrt(.Machine$double.eps))){
    stop("* PreEst.banded1 : 'upperK' should be an integer in [1,ncol(X)].")
  }
  #   3. delta :  Default 10, has to be larger than 2
  Banerjee14.delta = as.double(delta)
  if ((length(Banerjee14.delta)!=1)||(Banerjee14.delta<=2)||(is.na(Banerjee14.delta))||(is.infinite(Banerjee14.delta))){
    stop("* PreEst.banded1 : 'delta' should be a value larger than 2.")
  }
  #   4. logpi : log of prior distribution for bandwidth k.
  Banerjee14.logpi = logpi
  if (!is.function(Banerjee14.logpi)){
    stop("* PreEst.banded1 : 'logpi' should be a function.")
  }
  #   5. loss : type of loss used
  Banerjee14.loss = match.arg(loss)

  #-----------------------------------------------------
  ## MAIN COMPUTATION
  output = preest.Banerjee14(Banerjee14.upperK, X, Banerjee14.delta, Banerjee14.logpi, Banerjee14.loss)

  #-----------------------------------------------------
  ## RETURN
  return(output)
}
