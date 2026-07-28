# Original name : Fan13
#' Covariance Estimation via Hard Thresholding under Positive-Definiteness Constraint
#'
#' Sparse covariance estimation does not necessarily guarantee positive definiteness of an estimated
#' covariance matrix. Fan et al. (2013) proposed to solve this issue by taking an iterative procedure to
#' take an incremental decrease of threshold value until positive definiteness is preserved.
#'
#' The selected estimate is verified to have minimum eigenvalue at least
#' \eqn{\sqrt{\epsilon}}, where \eqn{\epsilon} is machine precision. Since
#' hard thresholding preserves the diagonal, the function stops with an
#' informative error when a column variance is below this tolerance and no
#' qualifying estimate can be constructed.
#'
#' @param X an \eqn{(n\times p)} matrix where each row is an observation.
#'
#' @return a named list containing: \describe{
#' \item{S}{a \eqn{(p\times p)} covariance matrix estimate.}
#' \item{optC}{a nonnegative threshold value \eqn{C_{min}} that guarantees
#' positive definiteness after thresholding. A value of zero means that the
#' sample covariance already satisfies the tolerance.}
#' }
#'
#' @examples
#' ## generate data from multivariate normal with Identity covariance.
#' pdim <- 5
#' data <- matrix(rnorm(10*pdim), ncol=pdim)
#'
#' ## apply 4 different schemes
#' out1 <- CovEst.hard(data, thr=0.1)  # threshold value 0.1
#' out2 <- CovEst.hard(data, thr=1)    # threshold value 1
#' out3 <- CovEst.hard(data, thr=10)   # threshold value 10
#' out4 <- CovEst.hardPD(data) # automatic threshold checking
#'
#' ## visualize 4 estimated matrices
#' mmessage <- paste("hardPD::optimal thr=",sprintf("%.2f",out4$optC),sep="")
#' gcol     <- gray((0:100)/100)
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
#' image(out1$S[,pdim:1], col=gcol, main="thr=0.1")
#' image(out2$S[,pdim:1], col=gcol, main="thr=1")
#' image(out3$S[,pdim:1], col=gcol, main="thr=10")
#' image(out4$S[,pdim:1], col=gcol, main=mmessage)
#' par(opar)
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
  X        = as.matrix(X)

  if ((!is.numeric(X))||(nrow(X)<2)||(ncol(X)<1)){
    stop("* CovEst.hardPD : X must be a numeric matrix with at least two rows and one column.")
  }

  sampleS = tryCatch(cov(X), error=function(e)e, warning=function(w)w)
  if (inherits(sampleS, "condition")||(!all(is.finite(sampleS)))){
    stop("* CovEst.hardPD : sample covariance computation failed.")
  }
  sampleS = (sampleS+t(sampleS))/2

  pdtol = sqrt(.Machine$double.eps)
  if (any(diag(sampleS)<pdtol)){
    stop("* CovEst.hardPD : every column must have variance at least sqrt(.Machine$double.eps).")
  }

  is_admissible = function(S){
    eigvals = tryCatch(eigen(S, symmetric=TRUE, only.values=TRUE)$values,
                       error=function(e)NULL)
    (!is.null(eigvals))&&all(is.finite(eigvals))&&(min(eigvals)>=pdtol)
  }
  threshold_covariance = function(thr){
    outS = sampleS
    offdiag = (row(outS)!=col(outS))
    outS[offdiag & (abs(outS)<=thr)] = 0
    return((outS+t(outS))/2)
  }

  #-----------------------------------------------------
  ## MAIN COMPUTATION
  # No thresholding is needed when the sample covariance already satisfies
  # the same positive-definiteness tolerance used by the search.
  if (is_admissible(sampleS)){
    message("* CovEst.hardPD : sample covariance itself is positive definite.")
    message("*               : So, we simply return Sample Covariance matrix.")
    return(list(S=sampleS, optC=0))
  }

  # Thresholding every off-diagonal entry gives a guaranteed fallback
  # whenever all marginal variances satisfy the tolerance.
  offdiag_values = sampleS[row(sampleS)!=col(sampleS)]
  Cmax = if (length(offdiag_values)>0){
    as.double(max(abs(offdiag_values)))
  } else {
    0
  }
  fallbackS = threshold_covariance(Cmax)
  if (!is_admissible(fallbackS)){
    stop("* CovEst.hardPD : no admissible positive-definite thresholded estimate exists.")
  }

  # Coarse scan from a known admissible threshold toward the unthresholded
  # covariance, retaining the last verified admissible candidate.
  coarse_grid = seq(from=Cmax, to=0, length.out=20)
  admissible_thr = Cmax
  admissible_S   = fallbackS
  inadmissible_thr = NA_real_
  for (i in 2:length(coarse_grid)){
    Ctmp = coarse_grid[i]
    tmpS = threshold_covariance(Ctmp)
    if (is_admissible(tmpS)){
      admissible_thr = Ctmp
      admissible_S   = tmpS
    } else {
      inadmissible_thr = Ctmp
      break
    }
  }

  # Fine scan moves from the known inadmissible endpoint toward the known
  # admissible endpoint and selects the first candidate that is admissible.
  Copt    = admissible_thr
  outputS = admissible_S
  if (is.finite(inadmissible_thr)){
    fine_grid = seq(from=inadmissible_thr, to=admissible_thr, length.out=20)
    for (Ctmp in fine_grid){
      tmpS = threshold_covariance(Ctmp)
      if (is_admissible(tmpS)){
        Copt    = Ctmp
        outputS = tmpS
        break
      }
    }
  }

  #-----------------------------------------------------
  ## RETURN OUTPUT
  # Guard against any unexpected numerical failure in the selected matrix.
  if (!is_admissible(outputS)){
    Copt    = Cmax
    outputS = fallbackS
  }
  if (!is_admissible(outputS)){
    stop("* CovEst.hardPD : failed to construct a positive-definite estimate.")
  }
  return(list(S=outputS, optC=as.double(Copt)))
}
