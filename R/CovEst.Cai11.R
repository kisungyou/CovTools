#' Covariance Estimation via Adaptive Thresholding
#'
#' Cai and Liu (2011) proposed an adaptive variant of Bickel and Levina (2008) - \code{\link{CovEst.Bickel08}}. The idea of \emph{adaptive thresholding} is
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
#' data <- mvtnorm::rmvnorm(100, sigma=diag(10))
#'
#' ## apply 4 different schemes
#' #  mthr is a vector of regularization parameters to be tested
#' mthr <- seq(from=0.01,to=0.99,length.out=10)
#'
#' out1 <- CovEst.Cai11(data, thr=0.1)  # threshold value 0.1
#' out2 <- CovEst.Cai11(data, thr=0.5)  # threshold value 0.5
#' out3 <- CovEst.Cai11(data, thr=0.5)  # threshold value 0.9
#' out4 <- CovEst.Cai11(data, thr=mthr) # automatic threshold checking
#'
#' ## visualize 4 estimated matrices
#' par(mfrow=c(2,2), pty="s")
#' image(pracma::flipud(out1$S), col=gray((0:100)/100), main="thr=0.1")
#' image(pracma::flipud(out2$S), col=gray((0:100)/100), main="thr=0.5")
#' image(pracma::flipud(out3$S), col=gray((0:100)/100), main="thr=0.9")
#' image(pracma::flipud(out4$S), col=gray((0:100)/100), main="automatic")
#'
#' @references
#' \insertRef{cai_adaptive_2011}{CovTools}
#'
#' @rdname CovEst_Cai11
#' @export
CovEst.Cai11 <- function(X, thr=0.5, nCV=10, parallel=TRUE){
  #-----------------------------------------------------
  ## PREPROCESSING
  fname    = "Cai11"
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
    output = covest.Cai11.once(X,thr)
  } else {
    output = covest.Cai11(X,nCV,nCore,thr)
  }

  #-----------------------------------------------------
  ## RETURN OUTPUT
  return(output)
}














#  ------------------------------------------------------------------------
#  Auxiliary Functions for CovEst.Cai11
#  ------------------------------------------------------------------------
# [Cai11] Adaptive Thresholding for Sparse CovMatEst ----------------------
#' @keywords internal
#' @noRd
covest.Cai11.once <- function(X,thrlevel){
  C = cor(X)
  C[which(abs(C)<=thrlevel)] = 0
  diag(C) = diag(cor(X))
  S = cov(X)
  dS2 = diag(sqrt(diag(S)))
  output = dS2 %*% C %*% dS2

  outputs = list()
  outputs$S = output
  # return outputs
  return(outputs)
}
#' @keywords internal
#' @noRd
covest.Cai11 <- function(X,nCV,nCore,thrvec){
  # get parameters
  n = nrow(X); p=ncol(X);
  nsearch = length(thrvec);
  # nCV
  S1 = array(0,c(p,p,nCV))
  S2 = array(0,c(p,p,nCV))
  nselect = round(n*(1-(1/log(n))))
  # divide and compute S for nfold cross validation
  idxall = (1:n)
  for (i in 1:nCV){
    idx1 = sample(idxall, nselect)
    idx2 = setdiff(idxall, idx1)
    S1[,,i] = cor(X[idx1,])
    S2[,,i] = cor(X[idx2,])
  }
  # total of 100 threshold values will be tested.
  stest  = sort(thrvec, decreasing=FALSE)
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
  C = cor(X)
  C[which(abs(C)<=thropt)] = 0
  diag(C) = diag(cor(X))
  S = cov(X)
  dS2 = diag(sqrt(diag(S)))
  output = dS2 %*% C %*% dS2


  outputs = list()
  outputs$S = output
  outputs$CV = data.frame(thr=stest, CVscore=as.vector(Rs))
  # return outputs
  return(outputs)
}

