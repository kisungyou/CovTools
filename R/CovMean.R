#' Estimate mean covariance matrix
#'
#' For a given 3-dimensional array where symmetric positive definite (SPD) matrices are stacked slice
#' by slice, it estimates Frechet mean on an open cone of SPD matrices under corresponding metric/distance
#' measure.
#'
#' @param A a \eqn{(p\times p\times N)} 3d array of \eqn{N} SPD matrices.
#' @param method the type of distance measures to be used; \code{"AIRM"} for Affine Invariant
#' Riemannian Metric,
#' \code{"Cholesky"} for Cholesky difference in Frobenius norm,
#' \code{"Euclidean"} for naive Frobenius norm as distance,
#' \code{"LERM"} for Log Euclidean Riemannian Metric,
#' \code{"Procrustes.SS"} for Procrustes Size and Shape measure,
#' \code{"Procrustes.Full"} for Procrustes analysis with scale,
#' \code{"PowerEuclidean"} for weighted eigenvalues by some exponent, and
#' \code{"RootEuclidean"} for matrix square root.
#' @param power a non-zero number for PowerEuclidean distance.
#' @return a \eqn{(p\times p)} mean covariance matrix estimated.
#'
#' @examples
#' ## generate 100 sample covariances of size (5-by-5).
#' samples = samplecovs(100,5)
#'
#' ## Compute mean of first 50 sample covariances from data under Normal(0,Identity).
#' mLERM = CovMean(samples[,,1:50],method="LERM")
#' mAIRM = CovMean(samples[,,1:50],method="AIRM")
#'
#' @references Dryden, I.L., Koloydenko, A., and Zhou, D. (2009) \emph{Non-Euclidean
#' statistics for covariance matrices, with applications to diffusion tensor imaging.} The
#' Annals of Applied Statistics, Vol.3(3):1102-1123.
#' @export
CovMean <- function(A,method=c("AIRM","Cholesky","Euclidean","LERM",
                               "Procrustes.SS","Procrustes.Full","PowerEuclidean",
                               "RootEuclidean"),power=1.0){
  ## PREPROCESSING
  ## 1) 3d array, 2) square, 3) symmetric, 4) sequentially check PDness
  if (length(dim(A))!=3){
    stop("* CovMean : 3d array should be used as input.")
  }
  if (dim(A)[1]!=dim(A)[2]){
    stop("* CovMean : A should be stack of square matrices.")
  }
  SymApply = apply(A, 3, isSymmetric, tol=sqrt(.Machine$double.eps))
  if (any(SymApply)==FALSE){
    SymIdx = which(!SymApply)
    if (length(SymIdx)==1){
      stopmessage = paste("* CovMean : slice number",SymApply,"is not Symmetric.")
    } else {
      stopmessage = paste("* CovMean : multiple of slices are not Symmetric.")
    }
    stop(stopmessage)
  }
  if (any(apply(A, 3, isSymmetric, tol=sqrt(.Machine$double.eps))==FALSE)){
    stop("* CovMean : each slice of A should be symmetric matrix.")
  }
  for (i in 1:dim(A)[3]){
    if (!check_pd(A[,,i])){
      stopmessage = paste(" CovMean : slice number",i,"is not Positive Definite.")
      stop(stopmessage)
    }
  }

  ## Main Iteration with Switch Argument
  if (missing(method)){method = "AIRM"}
  method = match.arg(method)
  if (method=="PowerEuclidean"){
    if (power==0){
      stop("* CovMean : 'power' should be a nonzero element. Suggests > 0.")
    }
    power = as.double(power)
  }
  outmean= switch(method,
                  Euclidean = meancov.Euclidean(A),
                  LERM = meancov.LogEuclidean(A),
                  Cholesky = meancov.Cholesky(A),
                  AIRM = meancov.Riemannian(A),
                  RootEuclidean=meancov.RootEuclidean(A),
                  Procrustes.SS=meancov.Procrustes.SS(A),
                  Procrustes.Full=meancov.Procrustes.Full(A),
                  PowerEuclidean= meancov.PowerEuclidean(A,power))
  return(outmean)
}
