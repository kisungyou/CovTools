#' Oracle Approximating Shrinkage Estimator
#'
#' Authors propose to estimate covariance matrix by iteratively approximating the shrinkage with
#' \deqn{\hat{\Sigma} = \rho \hat{F} + (1-\rho) \hat{S}}
#' where \eqn{\rho \in (0,1)} a control parameter/weight, \eqn{\hat{S}} an empirical covariance matrix, and \eqn{\hat{F}} a \emph{target} matrix.
#' It is proposed to use a structured estimate \eqn{\hat{F} = \textrm{Tr} (\hat{S}/p) \cdot I_{p\times p}} where \eqn{I_{p\times p}} is an identity matrix of dimension \eqn{p}.
#'
#' @param X an \eqn{(n\times p)} matrix where each row is an observation.
#'
#' @return a named list containing: \describe{
#' \item{S}{a \eqn{(p\times p)} covariance matrix estimate.}
#' \item{rho}{an estimate for convex combination weight.}
#' }
#'
#' @examples
#' ## CRAN-purpose small computation
#' # set a seed for reproducibility
#' set.seed(11)
#'
#' #  small data with identity covariance
#' pdim      <- 5
#' dat.small <- matrix(rnorm(10*pdim), ncol=pdim)
#'
#' #  run the code
#' out.small <- CovEst.2010OAS(dat.small)
#'
#' #  visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' image(diag(pdim)[,pdim:1],     main="true cov")
#' image(cov(dat.small)[,pdim:1], main="sample cov")
#' image(out.small$S[,pdim:1],    main="estimated cov")
#' par(opar)
#'
#' \dontrun{
#' ## want to see how delta is determined according to
#' #  the number of observations we have.
#' nsamples = seq(from=5, to=200, by=5)
#' nnsample = length(nsamples)
#'
#' #  we will record two values; rho and norm difference
#' vec.rho   = rep(0, nnsample)
#' vec.normd = rep(0, nnsample)
#' for (i in 1:nnsample){
#'   dat.norun <- matrix(rnorm(nsamples[i]*pdim), ncol=pdim) # sample in R^5
#'   out.norun <- CovEst.2010OAS(dat.norun)                  # run with default
#'
#'   vec.rho[i]   = out.norun$rho
#'   vec.normd[i] = norm(out.norun$S - diag(pdim),"f")       # Frobenius norm
#' }
#'
#' # let's visualize the results
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(nsamples, vec.rho,   lwd=2, type="b", col="red", main="estimated rhos")
#' plot(nsamples, vec.normd, lwd=2, type="b", col="blue",main="Frobenius error")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{chen_shrinkage_2010}{CovTools}
#'
#' @export
CovEst.2010OAS <- function(X){
  #-----------------------------------------------------
  ## PREPROCESSING
  fname    = "CovEst.2010RBLW"
  checker1 = invisible_datamatrix(X, fname)
  n = nrow(X) # number of observations
  p = ncol(X)

  #-----------------------------------------------------
  ## COMPUTATION
  #   1. MLE for Sigma
  Shat = stats::cov(X)*(nrow(X)-1)/(nrow(X))
  #      and some related values
  trS  = aux_trace(Shat)
  trS2 = aux_trace(Shat%*%Shat)

  #   2. structured estimate
  Fhat = (trS/p)*diag(p)

  #   3. rhohat
  term1 = (1-(2/p))*trS2 + (trS^2)        # IEEE publication is incorrect..
  term2 = (n+1-(2/p))*(trS2 - (trS^2)/p)  # but ArXiv version is correct...
  rhohat = max(min(term1/term2, 1), 0)

  #-----------------------------------------------------
  ## RETURN THE OUTPUT
  output = list()
  output$S   = (1-rhohat)*Shat + rhohat*Fhat
  output$rho = rhohat
  return(output)
}
