#' Covariance Estimation with Linear Shrinkage
#'
#' Ledoit and Wolf (2003, 2004) proposed a linear shrinkage strategy to estimate covariance matrix
#' with an application to portfolio optimization. An optimal covariance is written as a convex combination as follows,
#' \deqn{\hat{\Sigma} = \delta \hat{F} + (1-\delta) \hat{S}}
#' where \eqn{\delta \in (0,1)} a control parameter/weight, \eqn{\hat{S}} an empirical covariance matrix, and \eqn{\hat{F}} a \emph{target} matrix.
#' Although authors used \eqn{F} a highly structured estimator, we also enabled an arbitrary target matrix to be used as long as it's symmetric
#' and positive definite of corresponding size.
#'
#' @param X an \eqn{(n\times p)} matrix where each row is an observation.
#' @param target target matrix \eqn{F}. If \code{target=NULL}, \emph{constant correlation model} estimator is used. If \code{target} is specified as a qualified matrix, it is used instead.
#'
#' @return a named list containing: \describe{
#' \item{S}{a \eqn{(p\times p)} covariance matrix estimate.}
#' \item{delta}{an estimate for convex combination weight according to the relevant theory.}
#' }
#'
#' @examples
#' ## CRAN-purpose small computation
#' # set a seed for reproducibility
#' set.seed(11)
#'
#' #  small data with identity covariance
#' pdim      <- 5
#' dat.small <- matrix(rnorm(20*pdim), ncol=pdim)
#'
#' #  run the code with highly structured estimator
#' out.small <- CovEst.2003LW(dat.small)
#'
#' #  visualize
#' opar <- par(mfrow=c(1,3), pty="s")
#' image(diag(5)[,pdim:1], main="true cov")
#' image(cov(dat.small)[,pdim:1], main="sample cov")
#' image(out.small$S[,pdim:1], main="estimated cov")
#' par(opar)
#'
#' \dontrun{
#' ## want to see how delta is determined according to
#' #  the number of observations we have.
#' nsamples = seq(from=5, to=200, by=5)
#' nnsample = length(nsamples)
#'
#' #  we will record two values; delta and norm difference
#' vec.delta = rep(0, nnsample)
#' vec.normd = rep(0, nnsample)
#' for (i in 1:nnsample){
#'   dat.norun <- matrix(rnorm(nsamples[i]*pdim), ncol=pdim) # sample in R^5
#'   out.norun <- CovEst.2003LW(dat.norun)                   # run with default
#'
#'   vec.delta[i] = out.norun$delta
#'   vec.normd[i] = norm(out.norun$S - diag(pdim),"f")       # Frobenius norm
#' }
#'
#' # let's visualize the results
#' par(mfrow=c(1,2))
#' plot(nsamples, vec.delta, lwd=2, type="b", col="red", main="estimated deltas")
#' plot(nsamples, vec.normd, lwd=2, type="b", col="blue",main="Frobenius error")
#' }
#'
#' @references
#' \insertRef{ledoit_improved_2003}{CovTools}
#'
#' \insertRef{ledoit_well-conditioned_2004}{CovTools}
#'
#' \insertRef{ledoit_honey_2004}{CovTools}
#'
#' @rdname CovEst.2003LW
#' @export
CovEst.2003LW <- function(X, target=NULL){
  #-----------------------------------------------------
  ## PREPROCESSING
  # 1.
  fname    = "CovEst.2003LW"
  checker1 = invisible_datamatrix(X, fname)
  matS = stats::cov(X)*(nrow(X)-1)/(nrow(X))
  Y = t(X) # let's follow paper's orientation; (N x T) size
  parN = nrow(Y)
  parT = ncol(Y)

  # 2. target
  if ((length(target)==0)&&(is.null(target))){ # highly structured estimator
    matF = CovEst.2003LW.estF(matS)
  } else {                                     # user specified target matrix
    matF = target
  }

  #-----------------------------------------------------
  ## COMPUTATION : compute three parameters
  # 0. some common elements
  Ybar = as.vector(rowMeans(Y)) # (N) length
  matR = stats::cov2cor(matS)
  rbar = sum(matR[upper.tri(matR)])*2/(parN*(parN-1))
  # 1. pi hat
  matPi = covest_2003LW_computePi(Y,Ybar,matS)
  hatPi = sum(matPi)
  # 2. rho hat
  rho.term1 = sum(diag(matPi))
  rho.term2 = covest_2003LW_computeRho(Y, Ybar, matS, rbar)
  hatrho    = rho.term1 + rho.term2
  # 3. gamma hat
  hatgamma  = sum((matF-matS)^2)
  # 4. optimal delta value
  kappa = (hatPi - hatrho)/hatgamma
  delta = kappa/parT
  delta = max(min(delta, 1), 0)

  #-----------------------------------------------------
  ## COMPUTATION : squeezed matrix
  outS  = delta*matF + (1-delta)*matS

  #-----------------------------------------------------
  ## RETURN THE OUTPUT
  output = list()
  output$S = outS
  output$delta = delta
  return(output)
}



# auxiliary functions -----------------------------------------------------
#' @keywords internal
#' @noRd
CovEst.2003LW.estF <- function(matS){
  # parameters
  N = nrow(matS)

  # precompute
  matR = stats::cov2cor(matS)
  rbar = sum(matR[upper.tri(matR)])*2/(N*(N-1))
  matF = covest_2003LW_computeF(matS, rbar)
  return(matF)
}
