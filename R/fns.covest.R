# * Number of Parameters = 1 ----------------------------------------------
#  SINGLE CASE
#' @keywords internal
#' @noRd
thr1.once <- function(X, thr, func){
  output   = list()
  output$S = func(X, thr)
  return(output)
}
#' @keywords internal
#' @noRd
thr1.once.givenS <- function(S, thr, func_S){
  output   = list()
  output$S = func_S(S, thr)
  return(output)
}
#  MULTIPLE CASE
#' @keywords internal
#' @noRd
thr1.singlesum <- function(S1, S2, func_S, thr){
  N = dim(S1)[3]
  output = 0
  for (i in 1:N){
    S1tmp = func_S(S1[,,i], thr)
    S2tmp = S2[,,i]
    output = output + norm(S1tmp-S2tmp,"f")
  }
  return(output)
}
#' @keywords internal
#' @noRd
thr1.multiple <- function(X, nCV, nCore, func_S, thrvec){
  #-------------------------------------------------------------
  # 1. get parameters
  n = nrow(X)
  p = ncol(X)
  if ((length(nCV)!=1L)||(!is.numeric(nCV))||(!is.finite(nCV))||
      (nCV<1)||(nCV!=as.integer(nCV))){
    stop("* thr1.multiple : 'nCV' should be a positive integer.")
  }
  if ((length(nCore)!=1L)||(!is.numeric(nCore))||(!is.finite(nCore))||
      (nCore<1)||(nCore!=as.integer(nCore))){
    stop("* thr1.multiple : 'nCore' should be a positive integer.")
  }
  nCV  = as.integer(nCV)
  nCore = as.integer(nCore)
  # 2. nCV
  S1 = array(0,c(p,p,nCV))
  S2 = array(0,c(p,p,nCV))
  n1 = round(n*(1-(1/log(n))))
  n2 = round(n-n1)
  if ((!is.finite(n1))||(n1<2L)||(n2<2L)){
    stop(paste0("* thr1.multiple : the prescribed cross-validation split ",
                "requires at least two observations in both training and ",
                "validation samples."))
  }
  # 3. divide and compute S for nfold cross validation
  idxall = seq_len(n)
  for (i in seq_len(nCV)){
    idx1 = sample(idxall, n1)
    idx2 = setdiff(idxall, idx1)
    S1[,,i] = cov(X[idx1,,drop=FALSE])
    S2[,,i] = cov(X[idx2,,drop=FALSE])
  }
  # 4. covariance
  S = cov(X)

  #-------------------------------------------------------------
  # 1. get ready for multiple thresholding vectors
  stest   = thrvec
  nsearch = length(stest)
  if (nsearch<1L){
    stop("* thr1.multiple : 'thrvec' should contain at least one value.")
  }
  # 2. compute the scores sequentially unless parallel work was requested
  if (nCore==1){
    Rs = vapply(seq_len(nsearch), function(i){
      thr1.singlesum(S1, S2, func_S, stest[i])
    }, numeric(1))
  } else {
    cl = makeCluster(min(nCore, nsearch))
    on.exit(stopCluster(cl), add=TRUE)
    Rs = parallel::parSapply(
      cl,
      seq_len(nsearch),
      function(i, S1, S2, func_S, stest, score_fun){
        score_fun(S1, S2, func_S, stest[i])
      },
      S1=S1,
      S2=S2,
      func_S=func_S,
      stest=stest,
      score_fun=thr1.singlesum
    )
  }
  # 3. find the minimal element
  idxsmin = which(Rs==min(Rs))
  if (length(idxsmin)>1){
    idxsmin = idxsmin[1]
  }
  thropt = stest[idxsmin]

  #-------------------------------------------------------------
  # 1. compute the optimal one
  outputS = func_S(S, thropt)
  # 2. return output
  output = list()
  output$S  = outputS
  output$CV = data.frame(thr=stest, CVscore=as.vector(Rs))
  # 3. return output
  return(output)
}
