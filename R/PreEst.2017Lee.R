# Original name : Lee17
#' Bayesian Estimation of a Banded Precision Matrix (Lee 2017)
#'
#' \code{PreEst.2017Lee} returns a Bayes estimator of the banded precision matrix,
#' which is defined in subsection 3.3 of Lee and Lee (2017), using the k-BC prior.
#' The bandwidth is set at the mode of marginal posterior for the bandwidth parameter.
#'
#' @param X an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param upperK integer upper bound of bandwidth \eqn{k}, at most
#' \eqn{\min(p-1,n-3)}.
#' @param logpi log of prior distribution for bandwidth \eqn{k}. Default is a function proportional to \eqn{-k^4}.
#'
#' @return a named list containing: \describe{
#' \item{C}{a \eqn{(p\times p)} MAP estimate for precision matrix.}
#' }
#'
#' @examples
#' ## generate data from multivariate normal with Identity precision.
#' pdim = 5
#' data = matrix(rnorm(100*pdim), ncol=pdim)
#'
#' ## compare different K
#' out1 <- PreEst.2017Lee(data, upperK=1)
#' out2 <- PreEst.2017Lee(data, upperK=3)
#' out3 <- PreEst.2017Lee(data, upperK=4)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
#' image(diag(pdim)[,pdim:1], main="Original Precision")
#' image(out1$C[,pdim:1],     main="banded2::upperK=1")
#' image(out2$C[,pdim:1],     main="banded2::upperK=3")
#' image(out3$C[,pdim:1],     main="banded2::upperK=4")
#' par(opar)
#'
#' @references
#' \insertRef{lee_estimating_2017}{CovTools}
#'
#' @rdname PreEst.2017Lee
#' @export
PreEst.2017Lee <- function(X, upperK=min(floor(ncol(X)/2), ncol(X)-1, nrow(X)-3),
                           logpi=function(k){-k^4}){
  #-----------------------------------------------------
  ## PREPROCESSING
  #   1. typecheck : data
  fname    = "PreEst.2017Lee"
  checker1 = invisible_datamatrix(X, fname)
  if (!is.numeric(X)){
    stop("* PreEst.2017Lee : 'X' should be a numeric matrix.")
  }
  X = as.matrix(X)
  n = nrow(X)
  p = ncol(X)
  if ((n < 4)||(p < 2)){
    stop("* PreEst.2017Lee : 'X' should have at least 4 rows and 2 columns.")
  }
  #   2. upperK : upper bound of bandwidth
  maxK = min(p-1, n-3)
  if ((!is.numeric(upperK))||(length(upperK)!=1)||(!is.finite(upperK))||
      (upperK != floor(upperK))||(upperK<1)||(upperK>maxK)){
    stop(paste0("* PreEst.2017Lee : 'upperK' should be an integer in [1,",
                maxK, "]."))
  }
  Lee17.upperK = as.integer(upperK)
  #   3. logpi : log of prior distribution for bandwidth k.
  Lee17.logpi = logpi
  if (!is.function(Lee17.logpi)){
    stop("* PreEst.2017Lee : 'logpi' should be a function.")
  }

  #-----------------------------------------------------
  ## PREPROCESSING
  output      = preest.Lee17(Lee17.upperK, X, logpi=Lee17.logpi)

  #-----------------------------------------------------
  ## RETURN
  return(output)
}


# auxiliary functions -----------------------------------------------------
# [Lee.17] approximated posterior mean of banded precision matrix  --------
#######################################################################
# approximated posterior mean of banded precision matrix based on Chol. decomp. (Lee and Lee, 2017+)
# K : upper bound of bandwidth
# X : n-by-p data matrix
# logpi : log of prior distribution for bandwidth k. Default is a function proportional to -k^4.
#
# refer : Estimating Large Precision Matrices via Modified Cholesky Decomposition (https://arxiv.org/abs/1707.01143#)
#######################################################################
# main function
#' @keywords internal
#' @noRd
preest.Lee17<- function(K, X, logpi=function(k){ -k^4 }){
  X = scale(X, center = TRUE, scale =FALSE)

  n = nrow(X); p = ncol(X)
  A = matrix(0, p, p)
  Dinv = diag(1, nrow=p)

  logprob = vapply(seq_len(K), function(k){
    Lee17.logpostpi(k, X, logpi)
  }, numeric(1))
  if (any(!is.finite(logprob))){
    stop("* PreEst.2017Lee : posterior bandwidth score is not finite.")
  }
  best = max(logprob)
  tolerance = 1000*.Machine$double.eps*max(1, max(abs(logprob)))
  khatLL = as.integer(which(logprob >= (best-tolerance))[1]) # scalar posterior mode

  dhat1 = as.numeric(Lee17.dhat(1, khatLL, X))
  if ((length(dhat1)!=1)||(!is.finite(dhat1))||(dhat1<=0)){
    stop("* PreEst.2017Lee : residual variance estimate is not positive.")
  }
  Dinv[1, 1] = ((n-2)/(n*dhat1))
  for(j in 2:p){
    nj = n - min(j-1, khatLL) - 2
    ind = c(max(1,j-khatLL):(j-1))

    Xj = X[, j]
    Zj = Lee17.Z(j, khatLL, X)
    Covhat = t(Zj)%*%matrix(Xj)/n
    dhatj = as.numeric(Lee17.dhat(j, khatLL, X))
    if ((length(dhatj)!=1)||(!is.finite(dhatj))||(dhatj<=0)){
      stop("* PreEst.2017Lee : residual variance estimate is not positive.")
    }

    A[j, ind] = Lee17.solve(Lee17.VhatZ(j, khatLL, X), Covhat)
    Dinv[j, j] = ((nj)/(n*dhatj))
  }
  OmegaLL = (diag(1, nrow=p) - t(A))%*%Dinv%*%(diag(1, nrow=p) - A)
  OmegaLL = Lee17.ensurePD(OmegaLL)

  output = list()
  output$C = OmegaLL
  return(output)
}
#######################################################################
# auxiliary functions
#######################################################################
# Z_j(k) : n X k matrix
#' @keywords internal
#' @noRd
Lee17.Z <- function(j, k, X){
  ind = c(max(1, j-k):(j-1))
  res = X[, ind]
  if(j == 2) res = matrix(res)
  return(res)
}
# \hat{Var}(Z_j(k)) : k X k matrix
#' @keywords internal
#' @noRd
Lee17.VhatZ <- function(j, k, X){
  Zj = Lee17.Z(j, k, X)
  res = t(Zj)%*%Zj/nrow(X)
  return(res)
}
# \hat{d}_{jk}
#' @keywords internal
#' @noRd
Lee17.dhat <- function(j, k, X){
  Xj = X[, j]
  if(j == 1){
    res = sum(Xj^2)/nrow(X)   # scalar
  }else{
    Zj = Lee17.Z(j, k, X)
    VhatX = sum(Xj^2)/nrow(X)   # scalar
    Covhat = t(Zj)%*%matrix(Xj)/nrow(X)   # k X 1 matrix
    res = VhatX - t(Covhat)%*%Lee17.solve(Lee17.VhatZ(j, k, X), Covhat)
  }
  return(res)
}
# log of posterior probability for k : \log(\pi(k | \bfX_n)) (Lee and Lee, 2016)
#' @keywords internal
#' @noRd
Lee17.logpostpi <- function(k, X, logpi){
  p = ncol(X); n = nrow(X)
  logprior = logpi(k)
  if ((!is.numeric(logprior))||(length(logprior)!=1)||(!is.finite(logprior))){
    stop("* PreEst.2017Lee : 'logpi' should return one finite numeric value.")
  }
  res = logprior
  for(j in 2:p){
    nj = n - min(j-1,k) - 2
    Vdet = determinant(n*Lee17.VhatZ(j, k, X)/(2*pi), logarithm=TRUE)
    dhatj = as.numeric(Lee17.dhat(j, k, X))
    if ((Vdet$sign<=0)||(!is.finite(as.numeric(Vdet$modulus)))||
        (length(dhatj)!=1)||(!is.finite(dhatj))||(dhatj<=0)){
      stop("* PreEst.2017Lee : a bandwidth regression is singular.")
    }
    res = res - 0.5*as.numeric(Vdet$modulus) +
      lgamma(0.5*nj) - 0.5*nj*log(0.5*n*dhatj)
  }
  return(as.numeric(res))
}

#' @keywords internal
#' @noRd
Lee17.solve <- function(A, B=NULL){
  output = tryCatch({
    if (is.null(B)){
      solve(A)
    } else {
      solve(A, B)
    }
  }, error=function(e){
    stop("* PreEst.2017Lee : a regression design is singular.", call.=FALSE)
  })
  if (any(!is.finite(output))){
    stop("* PreEst.2017Lee : a regression solve produced non-finite values.")
  }
  return(output)
}

#' @keywords internal
#' @noRd
Lee17.ensurePD <- function(Omega){
  Omega = (Omega+t(Omega))/2
  if (any(!is.finite(Omega))){
    stop("* PreEst.2017Lee : precision estimate contains non-finite values.")
  }
  okay = tryCatch({
    chol(Omega)
    TRUE
  }, error=function(e) FALSE)
  if (!okay){
    stop("* PreEst.2017Lee : precision estimate is not positive definite.")
  }
  return(Omega)
}
