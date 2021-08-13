#' Compute Pairwise Distance for Symmetric Positive-Definite Matrices
#'
#' For a given 3-dimensional array where symmetric positive definite (SPD) matrices are stacked slice
#' by slice, it computes pairwise distance using various popular measures. Some of measures
#' are \emph{metric} as they suffice 3 conditions in mathematical context; nonnegative definiteness,
#' symmetry, and triangle inequalities. Other non-metric measures represent \emph{dissimilarities} between
#' two SPD objects.
#'
#' @param A a \eqn{(p\times p\times N)} 3d array of \eqn{N} SPD matrices.
#' @param method the type of distance measures to be used; \code{"AIRM"} for Affine Invariant
#' Riemannian Metric, \code{"Bhattacharyya"} for Bhattacharyya distance based on normal model,
#' \code{"Cholesky"} for Cholesky difference in Frobenius norm,
#' \code{"Euclidean"} for naive Frobenius norm as distance,
#' \code{"Hellinger"} for Hellinger distance based on normal model,
#' \code{"JBLD"} for Jensen-Bregman Log Determinant Distance,
#' \code{"KLDM"} for symmetrized Kullback-Leibler Distance Measure,
#' \code{"LERM"} for Log Euclidean Riemannian Metric,
#' \code{"Procrustes.SS"} for Procrustes Size and Shape measure,
#' \code{"Procrustes.Full"} for Procrustes analysis with scale,
#' \code{"PowerEuclidean"} for weighted eigenvalues by some exponent, and
#' \code{"RootEuclidean"} for matrix square root.
#'
#' @param power a non-zero number for PowerEuclidean distance.
#'
#' @return an \eqn{(N\times N)} symmetric matrix of pairwise distances.
#'
#' @examples
#' ## generate 100 SPD matrices of size (5-by-5)
#' samples = samplecovs(100,5)
#'
#' ## get pairwise distance for "AIRM"
#' distAIRM = CovDist(samples, method="AIRM")
#'
#' ## dimension reduction using MDS
#' ss = cmdscale(distAIRM)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' plot(ss[,1],ss[,2],main="2d projection")
#' par(opar)
#'
#' @references
#' \insertRef{arsigny_log-euclidean_2006}{CovTools}
#'
#' \insertRef{dryden_non-euclidean_2009}{CovTools}
#'
#' @export
CovDist <- function(A,method=c("AIRM","Bhattacharyya","Cholesky",
                               "Euclidean","Hellinger","JBLD","KLDM","LERM",
                               "Procrustes.SS","Procrustes.Full","PowerEuclidean",
                               "RootEuclidean"),power=1.0){
  ## PREPROCESSING
  ## 1) 3d array, 2) square, 3) symmetric, 4) sequentially check PDness
  if (length(dim(A))!=3){
    stop("* CovDist : 3d array should be used as input.")
  }
  if (dim(A)[1]!=dim(A)[2]){
    stop("* CovDist : A should be stack of square matrices.")
  }
  SymApply = apply(A, 3, isSymmetric, tol=sqrt(.Machine$double.eps))
  if (any(SymApply)==FALSE){
    SymIdx = which(!SymApply)
    if (length(SymIdx)==1){
      stopmessage = paste("* CovDist : slice number",SymApply,"is not Symmetric.")
    } else {
      stopmessage = paste("* CovDist : multiple of slices are not Symmetric.")
    }
    stop(stopmessage)
  }
  if (any(apply(A, 3, isSymmetric, tol=sqrt(.Machine$double.eps))==FALSE)){
    stop("* CovDist : each slice of A should be symmetric matrix.")
  }
  for (i in 1:dim(A)[3]){
    if (!check_pd(A[,,i])){
      stopmessage = paste(" CovDist : slice number",i,"is not Positive Definite.")
      stop(stopmessage)
    }
  }

  ## Main Iteration with Switch Argument
  if (missing(method)){method = "AIRM"}
  method = match.arg(method)
  if (method=="PowerEuclidean"){
    if (power==0){
      stop("* CovDist : 'power' should be a nonzero element. Suggests > 0.")
    }
    power = as.double(power)
  }
  outdist= switch(method,
                  AIRM = measure.AIRM.3d(A),
                  Bhattacharyya = measure.Bhattacharyya.3d(A),
                  Hellinger = measure.Hellinger.3d(A),
                  LERM = measure.LERM.3d(A),
                  KLDM = measure.KLDM.3d(A),
                  JBLD = measure.JBLD.3d(A),
                  Euclidean = measure.Euclidean.3d(A),
                  Cholesky = measure.Choleksy.3d(A),
                  RootEuclidean = measure.RootEuclidean.3d(A),
                  Procrustes.SS = measure.Procrustes.SS.3d(A),
                  Procrustes.Full=measure.Procrustes.Full.3d(A),
                  PowerEuclidean =measure.PowerEuclidean.3d(A,power)
                  )
  return(outdist)
}
