#  Original name : Qi06
#' Covariance Estimation via Nearest Positive-Definite Matrix Projection
#'
#' Qi and Sun (2006) proposed an algorithm for computing the positive correlation matrix
#' with Positive Definiteness and transforming it back in order to estimate covariance matrix.
#' This algorithm does not depend on any parameters.
#'
#' @param X an \eqn{(n\times p)} matrix where each row is an observation.
#'
#' @return a named list containing: \describe{
#' \item{S}{a \eqn{(p\times p)} covariance matrix estimate.}
#' }
#'
#' @examples
#' ## generate data from multivariate normal with Identity covariance.
#' pdim <- 5
#' data <- matrix(rnorm(10*pdim), ncol=pdim)
#'
#' ## compare against sample covariance
#' out1 <- cov(data)
#' out2 <- CovEst.nearPD(data) # apply nearPD
#'
#' ## visualize 2 estimated matrices
#' gcol <- gray((0:100)/100)
#' opar <- par(mfrow=c(1,2), pty="s")
#' image(out1[,pdim:1],   col=gcol, main="sample covariance")
#' image(out2$S[,pdim:1], col=gcol, main="SPD Projection")
#' par(opar)
#'
#' @references
#' \insertRef{qi_quadratically_2006}{CovTools}
#'
#' @rdname CovEst.nearPD
#' @export
CovEst.nearPD <- function(X){
  #-----------------------------------------------------
  ## PREPROCESSING
  fname    = "nearPD"
  checker1 = invisible_datamatrix(X, fname)

  #-----------------------------------------------------
  ## MAIN COMPUTATION
  #   1. covariance and correlation
  S     = cov(X)
  diagS = diag(S)
  C     = diag(1/sqrt(diagS))%*%S%*%diag(1/sqrt(diagS))

  #   2. apply nearPD
  Cadj  = matrix(Matrix::nearPD(C, corr=TRUE)$mat, nrow=nrow(C))

  #   3. adjust again to output
  outS  = diag(sqrt(diagS))%*%Cadj%*%diag(sqrt(diagS))

  #-----------------------------------------------------
  ## RETURN OUTPUT
  output = list()
  output$S = outS
  return(output)
}
