# Original name : Fan13
#' Covariance Estimation via Hard Thresholding under Positive-Definiteness Constraint
#'
#' Sparse covariance estimation does not necessarily guarantee positive definiteness of an estimated
#' covariance matrix. Fan et al. (2013) proposed to solve this issue by taking an iterative procedure to
#' take an incremental decrease of threshold value until positive definiteness is preserved.
#'
#' @param X an \eqn{(n\times p)} matrix where each row is an observation.
#'
#' @return a named list containing: \describe{
#' \item{S}{a \eqn{(p\times p)} covariance matrix estimate.}
#' \item{optC}{an optimal threshold value \eqn{C_{min}} that guarantees positive definiteness after thresholding.}
#' }
#'
#' @examples
#' ## generate data from multivariate normal with Identity covariance.
#' data <- mvtnorm::rmvnorm(5, sigma=diag(10))
#'
#' ## apply 4 different schemes
#' out1 <- CovEst.hard(data, thr=0.1)  # threshold value 0.1
#' out2 <- CovEst.hard(data, thr=1)    # threshold value 1
#' out3 <- CovEst.hard(data, thr=10)   # threshold value 10
#' out4 <- CovEst.hardPD(data) # automatic threshold checking
#'
#' ## visualize 4 estimated matrices
#' mmessage <- paste("hardPD::optimal thr=",sprintf("%.2f",out4$optC),sep="")
#' par(mfrow=c(2,2), pty="s")
#' image(pracma::flipud(out1$S), col=gray((0:100)/100), main="thr=0.1")
#' image(pracma::flipud(out2$S), col=gray((0:100)/100), main="thr=1")
#' image(pracma::flipud(out3$S), col=gray((0:100)/100), main="thr=10")
#' image(pracma::flipud(out4$S), col=gray((0:100)/100), main=mmessage)
#'
#' @references
#' \insertRef{fan_large_2013}{CovTools}
#'
#' @rdname CovEst.hardPD
#' @export
CovEst.hardPD <- function(X){
  #-----------------------------------------------------
  ## PREPROCESSING
  fname    = "CovEst.hardPD"
  checker1 = invisible_datamatrix(X, fname)

  #-----------------------------------------------------
  ## MAIN COMPUTATION
  # 1. find intermediate range
  tmpS = cov(X)
  diag(tmpS) = 0.0
  Cmax  = as.double(max(abs(tmpS)))
  Cvec1 = seq(from=Cmax, to=1e-8, length.out=20)
  iter  = 0
  for (i in 2:20){
    Ctmp = Cvec1[i]
    tmpS = CovEst.hard(X, thr=Ctmp)$S
    if ((min(eigen(tmpS, only.values=TRUE)$values)) < sqrt(.Machine$double.eps)){
      iter = i
      break
    }
  }
  if (iter==0){
    message("* CovEst.hardPD : sample covariance itself is positive definite.")
    message("*               : So, we simply return Sample Covariance matrix.")
    output = list()
    output$S = cov(X)
    return(output)
  } else {
    # 2. second stage
    Cvec2 = seq(from=Cvec1[iter], to=Cvec1[iter-1], length.out=20)
    Copt  = 0
    for (i in 1:20){
      Ctmp = Cvec2[i]
      tmpS = CovEst.hard(X, thr=Ctmp)$S
      if ((min(eigen(tmpS, only.values=TRUE)$values)) < sqrt(.Machine$double.eps)){
        Copt = Ctmp
        break
      }
    }

    #-----------------------------------------------------
    ## RETURN OUTPUT
    output = list()
    output$S = CovEst.hard(X, thr=Copt, nCV=1)$S
    output$optC = Copt
    return(output)
  }
}
