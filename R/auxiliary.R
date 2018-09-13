# generation --------------------------------------------------------------
#' Generate Sample Covariances of 2 groups
#'
#' For visualization purpose, \code{samplecovs} generates a 3d array
#' of stacked sample covariances where - in 3rd dimension, the first half
#' are sample covariances of samples generated independently from
#' normal distribution with identity covariance, where the latter half
#' consists of samples covariances from dense random population covariance.
#'
#' @param ncopy the total number of sample covariances to be generated.
#' @param size dimension \eqn{p}.
#'
#' @return a \eqn{(p\times p\times ncopy)} array of strictly positive definite sample covariances.
#' @examples
#' ## generate total of 20 samples covariances of size 5-by-5.
#' samples <- samplecovs(20,5)
#'
#'
#'@export
samplecovs <- function(ncopy,size){
  ## Parameter
  if ((length(ncopy)!=1)||(is.na(ncopy))||(is.infinite(ncopy))||(ncopy<3)){
    stop("* samplecovs : 'ncopy' should be a positive integer >= 3.")
  }
  if ((length(size)!=1)||(is.na(size))||(is.infinite(size))||(size<2)){
    stop("* samplecovs : 'size' should be a positive integer >= 2.")
  }
  ncopy = round(ncopy)
  size  = round(size)


  ## two groups
  xxx = rbind(matrix(rnorm(2*size*size),nrow=2*size),
              matrix(rnorm(2*size*size),nrow=2*size)+5)
  size1 = round(ncopy/2)      # number of first groups
  size2 = round(ncopy-size1)  # number of second groups
  sigma1 = diag(size)
  sigma2 = cov(xxx)

  covout = array(0,c(size,size,ncopy))
  for (i in 1:ncopy){
    if (i <= size1){
      sampledat = mvtnorm::rmvnorm(2*size, sigma=sigma1)
      covout[,,i] = cov(sampledat)
    } else {
      sampledat = mvtnorm::rmvnorm(2*size, sigma=sigma2)
      covout[,,i] = cov(sampledat)
    }
  }
  return(covout)
}


# CHECKERS ----------------------------------------------------------------
# 1. check if the given matrix is square sized
#' @keywords internal
#' @noRd
check_sqmat <- function(A){
  cond1 = (is.matrix(A)||(inherits(A, "Matrix")))
  cond2 = ((nrow(A)==ncol(A)))
  if (cond1&&cond2){
    return(TRUE)
  } else{
    return(FALSE)
  }
}
# 2. check if given matrix A is positive definite
#' @keywords internal
#' @noRd
check_pd <- function(A){
  eigens = eigen(A, symmetric=TRUE, only.values=TRUE)
  if (any(eigens$values<=0)){
    return(FALSE)
  } else {
    return(TRUE)
  }
}
# 3. check if a function is metric
#' @keywords internal
#' @noRd
check_metric <- function(func,n=5,thr=1e-5){
  n  = round(n)
  n2 = round(n^2)

  set.seed(496)
  A = matrix(rnorm(n2),ncol=n); A = cov(A);
  B = matrix(rnorm(n2),ncol=n); B = cov(B);
  C = matrix(rnorm(n2),ncol=n); C = cov(C);

  cond1 = (func(A,A)<thr)
  cond2 = ((abs(func(A,B)-func(B,A)))<thr)
  cond3 = (func(A,C)<=(func(A,B)+func(B,C)))

  if ((cond1)&&(cond2)&&(cond3)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
# 4. check if the given input is valid datamatrix without NA and Infs
#' @keywords internal
#' @noRd
check_datamatrix <- function(A){
  cond1 = (is.matrix(A)||(inherits(A, "Matrix")))
  cond2 = (!any(is.infinite(A)))
  cond3 = (!any(is.na(A)))
  if (cond1&&cond2&&cond3){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Get Series --------------------------------------------------------------
# 1. Get Inverse Root matrix
#' @keywords internal
#' @noRd
get_invroot <- function(X){
  eigs = eigen(X)
  if (any(eigs$values < .Machine$double.eps*10)){
    stop("** The desired covariance 'Sigma0' is invalid.")
  }
  out = eigs$vectors %*% diag((eigs$values)^(-0.5)) %*% t(eigs$vectors)
  return(out)
}
