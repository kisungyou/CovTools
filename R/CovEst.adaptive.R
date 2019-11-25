# Original name : Cai11
#' Covariance Estimation via Adaptive Thresholding
#'
#' Cai and Liu (2011) proposed an adaptive variant of Bickel and Levina (2008) - \code{\link{CovEst.hard}}. The idea of \emph{adaptive thresholding} is
#' to apply thresholding technique on correlation matrix in that it becomes \emph{adaptive} in terms of each variable.
#'
#' @param X an \eqn{(n\times p)} matrix where each row is an observation.
#' @param thr user-defined threshold value. If it is a vector of regularization values, it automatically selects one that minimizes cross validation risk.
#' @param nCV the number of repetitions for 2-fold random cross validations for each threshold value.
#' @param parallel a logical; \code{TRUE} to use half of available cores, \code{FALSE} to do every computation sequentially.
#'
#' @return a named list containing: \describe{
#' \item{S}{a \eqn{(p\times p)} covariance matrix estimate.}
#' \item{CV}{a dataframe containing vector of tested threshold values(\code{thr}) and corresponding cross validation scores(\code{CVscore}).}
#' }
#'
#' @examples
#' ## generate data from multivariate normal with Identity covariance.
#' pdim <- 5
#' data <- matrix(rnorm(10*pdim), ncol=pdim)
#'
#' ## apply 4 different schemes
#' #  mthr is a vector of regularization parameters to be tested
#' mthr <- seq(from=0.01,to=0.99,length.out=10)
#'
#' out1 <- CovEst.adaptive(data, thr=0.1)  # threshold value 0.1
#' out2 <- CovEst.adaptive(data, thr=0.5)  # threshold value 0.5
#' out3 <- CovEst.adaptive(data, thr=0.5)  # threshold value 0.9
#' out4 <- CovEst.adaptive(data, thr=mthr) # automatic threshold checking
#'
#' ## visualize 4 estimated matrices
#' opar <- par(mfrow=c(2,2), pty="s")
#' image(out1$S[,pdim:1], col=gray((0:100)/100), main="thr=0.1")
#' image(out2$S[,pdim:1], col=gray((0:100)/100), main="thr=0.5")
#' image(out3$S[,pdim:1], col=gray((0:100)/100), main="thr=0.9")
#' image(out4$S[,pdim:1], col=gray((0:100)/100), main="automatic")
#' par(opar)
#'
#' @references
#' \insertRef{cai_adaptive_2011}{CovTools}
#'
#' @rdname CovEst.adaptive
#' @export
CovEst.adaptive <- function(X, thr=0.5, nCV=10, parallel=FALSE){
  #-----------------------------------------------------
  ## PREPROCESSING
  fname    = "CovEst.adaptive"
  pnameTHR = "'thr'"
  pnamenCV = "'nCV"
  pnamenthrs = "'nthrs'"
  #   1. data matrix
  checker1 = invisible_datamatrix(X, fname)
  #   2. thr
  if (length(as.vector(thr))==1){
    checker3 = invisible_PosReal(thr, fname, pnameTHR)
    CV = FALSE
  } else { # vector threshold value case
    nthrs = length(thr)
    for (i in 1:nthrs){
      checker3 = invisible_PosReal(thr[i], fname, pnameTHR)
    }
    CV = TRUE
  }
  #   4. nCV
  if (CV==TRUE){
    nCV = as.integer(nCV)
    checker4 = invisible_PosIntMM(nCV, fname, pnamenCV, 1, nrow(X))
  }
  #   5. parallel
  if (!parallel){
    nCore = 1
  } else {
    nCore = max(round(detectCores()/2),1)
  }

  #-----------------------------------------------------
  ## MAIN COMPUTATION
  if (CV==FALSE){
    output = thr1.once(X,thr,func_adaptive)
  } else {
    output = thr1.multiple(X,nCV,nCore,func_adaptive_givenS,thr)
  }

  #-----------------------------------------------------
  ## RETURN OUTPUT
  return(output)
}


# Auxiliary Functions for Adaptive Thresholding ---------------------------
#' @keywords internal
#' @noRd
func_adaptive <- function(X, thr){
  #   1. compute S and diagonal term
  S     = cov(X)
  diagS = diag(S)

  #   2. transform it into correlation matrix
  R = diag(1/sqrt(diagS))%*%S%*%diag(1/sqrt(diagS))

  #   3. apply bickel's hard thresholding
  outR  = R
  diagR = diag(R)
  outR[which(abs(outR)<=thr)] = 0
  diag(outR) = diagR

  #   4. return output
  outS = diag(sqrt(diagS))%*%outR%*%diag(sqrt(diagS))

  #   3. return output
  return(outS)
}
#' @keywords internal
#' @noRd
func_adaptive_givenS <- function(S, thr){
  #   1. compute S and diagonal term
  diagS = diag(S)

  #   2. transform it into correlation matrix
  R = diag(1/sqrt(diagS))%*%S%*%diag(1/sqrt(diagS))

  #   3. apply bickel's hard thresholding
  outR  = R
  diagR = diag(R)
  outR[which(abs(outR)<=thr)] = 0
  diag(outR) = diagR

  #   4. return output
  outS = diag(sqrt(diagS))%*%outR%*%diag(sqrt(diagS))

  #   3. return output
  return(outS)
}
