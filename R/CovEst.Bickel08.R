#' Covariance Estimation via Thresholding
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
#' data <- mvtnorm::rmvnorm(100, sigma=diag(5))
#'
#' ## apply 4 different schemes
#' #  mthr is a vector of regularization parameters to be tested
#' mthr <- exp(seq(from=log(0.1),to=log(10),length.out=10))
#'
#' out1 <- CovEst.Bickel08(data, thr=0.1)  # threshold value 0.1
#' out2 <- CovEst.Bickel08(data, thr=1)    # threshold value 1
#' out3 <- CovEst.Bickel08(data, thr=10)   # threshold value 10
#' out4 <- CovEst.Bickel08(data, thr=mthr) # automatic threshold checking
#'
#' ## visualize 4 estimated matrices
#' par(mfrow=c(2,2), pty="s")
#' image(pracma::flipud(out1$S), col=gray((0:100)/100), main="thr=0.1")
#' image(pracma::flipud(out2$S), col=gray((0:100)/100), main="thr=1")
#' image(pracma::flipud(out3$S), col=gray((0:100)/100), main="thr=10")
#' image(pracma::flipud(out4$S), col=gray((0:100)/100), main="automatic")
#'
#' @references
#' \insertRef{bickel_covariance_2008}{CovTools}
#'
#' @rdname CovEst_Bickel08
#' @export
CovEst.Bickel08 <- function(X, thr=sqrt(log(ncol(X))/nrow(X)), nCV=10, parallel=TRUE){
  #-----------------------------------------------------
  ## PREPROCESSING
  fname    = "Bickel08"
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
  checker4 = invisible_PosIntMM(nCV, fname, pnamenCV, 1, nrow(X))
  #   5. parallel
  if (!parallel){
    nCore = 1
  } else {
    nCore = max(round(detectCores()/2),1)
  }

  #-----------------------------------------------------
  ## MAIN COMPUTATION
  if (CV==FALSE){
    output = covest.Bickel08.once(X,thr)
  } else {
    output = covest.Bickel08(X,nCV,nCore,thr)
  }

  #-----------------------------------------------------
  ## RETURN OUTPUT
  return(output)
}





#  ------------------------------------------------------------------------
#  Auxiliary Functions for CovEst.Bickel08
#  ------------------------------------------------------------------------
# [Bickel08] Covariance Regularization by Thresholding --------------------
#' @keywords internal
#' @noRd
covest.Bickel08.once <- function(X,thrlevel){
  S = cov(X)
  outputS = S

  outputS[which(abs(S)<=thrlevel)] = 0
  diag(outputS) = diag(S)

  output = list()
  output$S = outputS
  return(output)
}
#' @keywords internal
#' @noRd
covest.Bickel08 <- function(X,nCV,nCore,thrvec){
  # get parameters
  n = nrow(X); p=ncol(X);
  # nCV
  S1 = array(0,c(p,p,nCV))
  S2 = array(0,c(p,p,nCV))
  n1 = round(n*(1-(1/log(n))))
  n2 = round(n-n1)
  # divide and compute S for nfold cross validation
  idxall = (1:n)
  for (i in 1:nCV){
    idx1 = sample(idxall, n1)
    idx2 = setdiff(idxall, idx1)
    S1[,,i] = cov(X[idx1,])
    S2[,,i] = cov(X[idx2,])
  }
  # Compute Sample Covariance anyway
  S = cov(X)
  # now it's changed
  #   we are given
  stest   = thrvec
  nsearch = length(stest)
  # parallel setup
  if (nCore==1){
    cl = makeCluster(1)
    registerDoParallel(cl)
  } else {
    cl = makeCluster(nCore)
    registerDoParallel(cl)
  }
  # parallel computation let's go
  # at first, let's not use
  itforeach=NULL
  Rs = foreach (itforeach=1:nsearch, .combine=cbind) %dopar% {
    covest.Bickel08.singlesum(S1,S2,stest[itforeach])
  }
  # stop cluster
  stopCluster(cl)
  # find the minimal element.
  idxsmin = which(Rs==min(Rs))
  if (length(idxsmin)>1){
    idxsmin = idxsmin[1]
  }
  # optimal threshold value
  thropt = stest[idxsmin]
  # Compute Sample Covariance anyway
  outputS = S
  outputS[which(abs(S)<=thropt)] = 0
  diag(outputS) = diag(S)

  output = list()
  output$S = outputS
  output$CV = data.frame(thr=stest, CVscore=as.vector(Rs))
  # return output
  return(output)
}
#' @keywords internal
#' @noRd
covest.Bickel08.singlesum <- function(S1, S2, thr){
  N = dim(S1)[3]
  output = 0
  for (i in 1:N){
    S1tmp = S1[,,i]
    S2tmp = S2[,,i]
    S1tmp[which(abs(S1tmp)<=thr)] = 0
    diag(S1tmp) = diag(S1[,,i])
    output = output + norm(S1tmp-S2tmp, "f")
  }
  return(output)
}
