# [2009.Dryden] Main 1. Euclidean -----------------------------------------
#' @keywords internal
#' @noRd
meancov.Euclidean <- function(A){
  output = apply(A, c(1,2), mean)
  return(output)
}

# [2009.Dryden] Main 2. LogEuclidean --------------------------------------
#' @keywords internal
#' @noRd
meancov.LogEuclidean <- function(A){
  p = dim(A)[1]
  N = dim(A)[3]
  output = array(0,c(p,p))
  # log sum
  for (i in 1:N){
    output = expm::logm(A[,,i])+output
  }
  output = expm::expm(output/N)
  return(output)
}

# [2009.Dryden] Main 3. Cholesky ------------------------------------------
#' @keywords internal
#' @noRd
meancov.Cholesky <- function(A){
  p = dim(A)[1]
  N = dim(A)[3]
  outChol = array(0,c(p,p,N))
  # 3-1. Cholesky Decomposition
  for (i in 1:N){
    tmp= tryCatch(t(chol(A[,,i])), error=function(e)e, warning=function(w)w)
    if (inherits(tmp, "error")){
      stop("* CovMean.Cholesky : cholesky decomposition failed.")
    } else if (inherits(tmp, "warning")){
      message("* CovMean.Cholesky : warning")
      message(paste("* CovMean.Cholesky :"),tmp$warning)
      return(NA)
    } else {
      outChol[,,i] = tmp
    }
  }
  # 3-2. averaging
  outDelta = apply(outChol, c(1,2), mean)
  # 3-3. output
  output = outDelta%*%t(outDelta)
  return(output)
}

# [2009.Dryden] Main 4. Riemannian (or, AIRM) -----------------------------
#' @keywords internal
#' @noRd
meancov.Riemannian <- function(A){
  output = port_estcov(A, method="AIRM")
  return(output)
}

# [2009.Dryden] Main 5. Root Euclidean ------------------------------------
#' @keywords internal
#' @noRd
meancov.RootEuclidean <- function(A){
  p = dim(A)[1]
  N = dim(A)[3]
  outRoot = array(0,c(p,p,N))
  # 5-1. square root of matrix
  for (i in 1:N){
    tmp = tryCatch(expm::sqrtm(A[,,i]), error=function(e)e, warning=function(w)w)
    if (inherits(tmp, "error")){
      stop("* CovMean.RootEuclidean : acquiring matrix square root failed.")
    } else if (inherits(tmp, "warning")){
      message("* CovMean.RootEuclidean : warning")
      message(paste("* CovMean.RootEuclidean :"),tmp$warning)
      return(NA)
    } else {
      outRoot[,,i] = tmp
    }
  }
  # 5-2. acquire Delta_H
  outDelta = apply(outRoot, c(1,2), mean)
  # 5-3. output
  output = outDelta%*%t(outDelta)
  return(output)
}


# [2009.Dryden] Main 6. Procrustes Size and Shape -------------------------
# NOTE : 'shapes' requires 'rgl', which requires following instructions.
# https://stackoverflow.com/questions/31820865/error-in-installing-rgl-package
# Also, use 'port_estcov'
#' @keywords internal
#' @noRd
meancov.Procrustes.SS <- function(A){
  output = port_estcov(A, method="Procrustes.SS")
  return(output)
}


# [2009.Dryden] Main 7. Procrustes Full -----------------------------------
#' @keywords internal
#' @noRd
meancov.Procrustes.Full <- function(A){
  output = port_estcov(A, method="Procrustes.Full")
  return(output)
}


# [2009.Dryden] Main 8. Power Euclidean -----------------------------------
#' @keywords internal
#' @noRd
meancov.PowerEuclidean <- function(A, pcoef){
  # 1.
  p = dim(A)[1]
  M = dim(A)[3]
  # 2. S_i^\alpha
  Salpha = array(0,dim(A))
  for (i in 1:M){
    eigtgt = eigen(A[,,i], symmetric=TRUE)
    Salpha[,,i] = eigtgt$vectors %*% diag((eigtgt$values)^pcoef) %*% t(eigtgt$vectors)
  }
  # 3. mean estimation
  outDelta = apply(Salpha, c(1,2), mean)
  eigDelta = eigen(outDelta)
  output   = eigDelta$vectors %*% (diag(eigDelta$values)^(1/pcoef)) %*% t(eigDelta$vectors)
  return(output)
}
