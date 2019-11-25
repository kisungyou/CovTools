# Original Name : Bickel08
#' Covariance Estimation via Hard Thresholding
#'
#' Bickel and Levina (2008) proposed a sparse covariance estimation technique to apply thresholding on off-diagonal elements of
#' the sample covariance matrix. The entry of sample covariance matrix \eqn{S_{i,j}=0} if \eqn{|S_{i,j}|<=\tau} where \eqn{\tau} is
#' a thresholding value (\code{thr}). If \code{thr} is rather a vector of regularization parameters, it applies
#' cross-validation scheme to select an optimal value.
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
#' mthr <- exp(seq(from=log(0.1),to=log(10),length.out=10))
#'
#' out1 <- CovEst.hard(data, thr=0.1)  # threshold value 0.1
#' out2 <- CovEst.hard(data, thr=1)    # threshold value 1
#' out3 <- CovEst.hard(data, thr=10)   # threshold value 10
#' out4 <- CovEst.hard(data, thr=mthr) # automatic threshold checking
#'
#' ## visualize 4 estimated matrices
#' gcol <- gray((0:100)/100)
#' opar <- par(mfrow=c(2,2), pty="s")
#' image(out1$S[,pdim:1], col=gcol, main="thr=0.1")
#' image(out2$S[,pdim:1], col=gcol, main="thr=1")
#' image(out3$S[,pdim:1], col=gcol, main="thr=10")
#' image(out4$S[,pdim:1], col=gcol, main="automatic")
#' par(opar)
#'
#' @references
#' \insertRef{bickel_covariance_2008}{CovTools}
#'
#' @rdname CovEst.hard
#' @export
CovEst.hard <- function(X, thr=sqrt(log(ncol(X))/nrow(X)), nCV=10, parallel=FALSE){
  #-----------------------------------------------------
  ## PREPROCESSING
  fname    = "CovEst.hard"
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
  nCV = as.integer(nCV)
  if (CV==TRUE){
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
    output = thr1.once(X,thr,func_hard)
  } else {
    output = thr1.multiple(X,nCV,nCore,func_hard_givenS,thr)
  }

  #-----------------------------------------------------
  ## RETURN OUTPUT
  return(output)
}


# Auxiliary Functions for Hard Thresholding -------------------------------
#' @keywords internal
#' @noRd
func_hard <- function(X, thr){
  #   1. compute S and diagonal term
  S     = cov(X)
  outS  = S
  diagS = diag(outS)

  #   2. adjust
  outS[which(abs(S)<=thr)] = 0
  diag(outS) = diagS

  #   3. return output
  return(outS)
}
#' @keywords internal
#' @noRd
func_hard_givenS <- function(S, thr){
  #   1. compute S and diagonal term
  outS  = S
  diagS = diag(outS)

  #   2. adjust
  outS[which(abs(S)<=thr)] = 0
  diag(outS) = diagS

  #   3. return output
  return(outS)
}
