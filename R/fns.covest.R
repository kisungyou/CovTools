


# [Donoho95] Wavelet Shrinkage: Asymptopia? -------------------------------
#' @keywords internal
#' @noRd
covest.Donoho95.once <- function(X,thrlevel){
  # Apply Soft Thresholding
  S    = cov(X)
  outS = covest.Donoho95.applysoft(S, thrlevel)
  output = list()
  output$S = outS
  return(output)
}
#' @keywords internal
#' @noRd
covest.Donoho95 <- function(X,nCV,nCore,nsearch){
  # get parameters
  n = nrow(X); p=ncol(X);
  # nCV
  S1 = array(0,c(p,p,nCV))
  S2 = array(0,c(p,p,nCV))
  nselect = round(n*(1-(1/log(n))))
  # divide and compute S for nfold cross validation
  idxall = (1:n)
  for (i in 1:nCV){
    idx1 = sample(idxall, nselect)
    idx2 = setdiff(idxall, idx1)
    S1[,,i] = cov(X[idx1,])
    S2[,,i] = cov(X[idx2,])
  }
  # total of 100 threshold values will be tested.
  S = cov(X)
  stest  = seq(from=(0.01*max(abs(S))),to=(0.99*max(abs(S))),length.out=nsearch)
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
    covest.Donoho95.singlesum(S1,S2,stest[itforeach])
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
  # Apply Soft Thresholding
  outS = covest.Donoho95.applysoft(S, thropt)
  output = list()
  output$S = outS
  output$CV = data.frame(thr=stest, CVscore=as.vector(Rs))
  return(output)
}
#' @keywords internal
#' @noRd
covest.Donoho95.singlesum <- function(S1, S2, thr){
  N = dim(S1)[3]
  output = 0
  for (i in 1:N){
    S1tmp = covest.Donoho95.applysoft(S1[,,i], thr)
    S2tmp = S2[,,i]
    output= output + norm(S1tmp-S2tmp, "f")
  }
  return(output)
}
#' @keywords internal
#' @noRd
covest.Donoho95.applysoft <- function(S, thr){
  output = array(0,dim(S))
  idxPos = which(S>0)
  idxNeg = which(S<0)

  # positive part
  output[idxPos] = pmax(S[idxPos]-thr,0)
  output[idxNeg] = pmin(0,S[idxNeg]+thr)
  diag(output)   = diag(S)
  return(output)
}


# [Fan.2013] Large covariance estimation by thresholding pc ---------------
#' @keywords internal
#' @noRd
covest.Fan13.once <- function(X,thropt){
  S    = cov(X)
  outS = covest.Fan13.singlethr(S, thropt)
  output = list()
  output$S = outS
  return(output)
}
#' @keywords internal
#' @noRd
covest.Fan13 <- function(X,nCV,nCore,nsearch){
  # get parameters
  n = nrow(X); p=ncol(X);
  # nCV
  S1 = array(0,c(p,p,nCV))
  S2 = array(0,c(p,p,nCV))
  nselect = round(n*(1-(1/log(n))))
  # divide and compute S for nfold cross validation
  idxall = (1:n)
  for (i in 1:nCV){
    idx1 = sample(idxall, nselect)
    idx2 = setdiff(idxall, idx1)
    S1[,,i] = cov(X[idx1,])
    S2[,,i] = cov(X[idx2,])
  }

  ## 1st. we need to find Cmin
  covX   = cov(X)
  SS     = covX
  Cflag  = TRUE
  Cmax = max(abs(setdiff(as.vector(SS),diag(covX))))
  Cmin = Cmax
  for (i in 1:1000){
    Cmin   = 0.99*Cmin  #update Cmax
    SS     = covest.Fan13.singlethr(SS, Cmin)
    eigSS  = min(eigen(SS, only.values = TRUE)$values)
    if (eigSS<=sqrt(sqrt(.Machine$double.eps))){
      Cmin = Cmin/0.99
      break
    }
  }
  if (Cmin==Cmax){
    message("* CovEst.Fan13 : Cmin is identical to Cmax.")
    Cmin = 0.9999*Cmax
  }

  ## 2nd. now, crossvalidation is run.
  stest  = seq(from=Cmin, to=Cmax, length.out=nsearch)
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
  # Apply Hard Thresholding
  S    = cov(X)
  outS = covest.Fan13.singlethr(S, thropt)
  output = list()
  output$S = outS
  output$CV = data.frame(thr=stest, CVscore=as.vector(Rs))
  return(output)
}
#' @keywords internal
#' @noRd
covest.Fan13.singlethr <- function(SS, thr){
  S = SS
  S[which(abs(S)<=thr)] = 0
  diag(S) = diag(SS)
  return(S)
}


# [Qi.2006] nearPD --------------------------------------------------------
#' @keywords internal
#' @noRd
covest.Qi06.once <- function(X){
  diagShalf = diag(sqrt(diag(cov(X))))
  Rt        = cor(X)
  Rhat      = matrix(Matrix::nearPD(Rt,corr=TRUE,keepDiag=TRUE)$mat, nrow=nrow(Rt))
  Shat      = diagShalf%*%Rhat%*%diagShalf

  output    = list()
  output$S  = Shat
  return(output)
}
#' @keywords internal
#' @noRd
covest.Qi06 <- function(X){
  diagShalf = diag(sqrt(diag(cov(X))))
  Rt        = cor(X)
  Rhat      = matrix(Matrix::nearPD(Rt,corr=TRUE,keepDiag=TRUE)$mat, nrow=nrow(Rt))
  Shat      = diagShalf%*%Rhat%*%diagShalf

  output    = list()
  output$S  = Shat
  output$CV = "* CovEst : nearest correlation matrix method does not involve cross validation."
  return(output)
}

