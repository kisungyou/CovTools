################################################################
## COLLECTION OF DISTANCE MEASURES IN 3D
################################################################
# 1. AIRM -----------------------------------------------------------------
#' Metric based on Generalized Eigenvalue Decomposition / Affine Invariant Riemannian Metric
#' @keywords internal
#' @noRd
measure.AIRM <- function(A,B){
  ## three conditions : size argument, symmetric, positive definite
  # 1. sqmat
  if ((!check_sqmat(A))||(!check_sqmat(B))){
    stop("* measure.AIRM : input is not square matrix/Matrix.")
  }
  # 2. symmetric
  leveltol = sqrt(.Machine$double.eps)
  if ((!isSymmetric(A, tol=leveltol))||(!isSymmetric(B, tol=leveltol))){
    stop("* measure.AIRM : input is not symmetric.")
  }
  A = (A+t(A))/2
  B = (B+t(B))/2
  # 3. positive definite
  if ((!check_pd(A))||(!check_pd(B))){
    stop("* measure.AIRM : input is not positive definite.")
  }

  ## Main Computation
  objeig = tryCatch(geigen(A, B, symmetric=TRUE, only.values = TRUE), error=function(e)e, warning=function(w)w)
  if (inherits(objeig, "error")){
    stop("* measure.AIRM : generalized eigenvalue decomposition failed.")
  } else if (inherits(objeig, "warning")){
    message(paste("* measure.AIRM :",objeig$message))
    return(NA)
  } else {
    output = sqrt(sum((log(objeig$values))^2))
    return(output)
  }
}
#' @keywords internal
#' @noRd
measure.AIRM.3d <- function(array3d){
  # 1. get size and set ready
  M      = dim(array3d)[3]
  data3d = array(0,dim(array3d))
  # 2. a priori symmetrization
  for (i in 1:M){
    data3d[,,i] = (array3d[,,i]+t(array3d[,,i]))/2
  }
  # 3. main computation using geigen
  outdist = array(0,c(M,M))
  for (i in 1:(M-1)){
    A = data3d[,,i]
    for (j in (i+1):M){
      B = data3d[,,j]
      objeig = tryCatch(geigen(A, B, symmetric=TRUE, only.values = TRUE), error=function(e)e, warning=function(w)w)
      if (inherits(objeig, "error")){
        stop("* measure.AIRM : generalized eigenvalue decomposition failed.")
      } else if (inherits(objeig, "warning")){
        message(paste("* measure.AIRM :",objeig$message))
        return(NA)
      } else {
        outdist[i,j] = sqrt(sum((log(objeig$values))^2))
        outdist[j,i] = sqrt(sum((log(objeig$values))^2))
      }
    }
  }
  return(outdist)
}

# 2. Bhattacharyya --------------------------------------------------------
#' Bhattacharyya distance with Normal Assumption : DISTANCE
#' @keywords internal
#' @noRd
measure.Bhattacharyya.3d <- function(array3d){
  # 1. get parameters
  M = dim(array3d)[3]
  # 2. determinant computation for each
  det3d = array(0,c(1,M))
  for (i in 1:M){
    detA = tryCatch(det(array3d[,,i]), error=function(e)e, warning=function(w)w)
    if (inherits(detA, "error")){
      stop("* measure.Bhattacharyya : computing determinant failed.")
    } else if (inherits(detA, "warning")){
      message(paste("* measure.Bhattacharyya :",detA$message))
      return(NA)
    } else {
      det3d[i] = detA
    }
  }
  # 3. main iteration
  outdist = array(0,c(M,M))
  for (i in 1:(M-1)){
    A    = array3d[,,i]
    detA = det3d[i]
    for (j in (i+1):M){
      B    = array3d[,,j]
      detB = det3d[j]
      C    = (A+B)/2
      detC = tryCatch(det(C), error=function(e)e, warning=function(w)w)
      if (inherits(detC, "error")){
        stop("* measure.Bhattacharyya : computing determinant failed.")
      } else if (inherits(detC, "warning")){
        message(paste("* measure.Bhattacharyya :",detA$message))
        return(NA)
      }
      value = log(detC/sqrt(detA*detB))/2
      outdist[i,j] = value
      outdist[j,i] = value
    }
  }
  return(outdist)
}

# 3. Hellinger ------------------------------------------------------------
#' Hellinger Distance on normal model : METRIC
#' @keywords internal
#' @noRd
measure.Hellinger.3d <- function(array3d){
  M       = dim(array3d)[3]
  DBmat   = measure.Bhattacharyya.3d(array3d)
  outdist = array(0, c(M,M))
  for (i in 1:(M-1)){
    for (j in (i+1):M){
      outdist[i,j] = sqrt(1-exp(-DBmat[i,j]))
      outdist[j,i] = sqrt(1-exp(-DBmat[i,j]))
    }
  }
  return(outdist)
}

# 4. LERM -----------------------------------------------------------------
#' Log-Euclidean Riemannian Metric
#' @keywords internal
#' @noRd
measure.LERM <- function(A,B){
  ## three conditions : size argument, symmetric, positive definite
  # 1. sqmat
  if ((!check_sqmat(A))||(!check_sqmat(B))){
    stop("* measure.LERM : input is not square matrix/Matrix.")
  }
  # 2. symmetric
  leveltol = sqrt(.Machine$double.eps)
  if ((!isSymmetric(A, tol=leveltol))||(!isSymmetric(B, tol=leveltol))){
    stop("* measure.LERM : input is not symmetric.")
  }
  A = (A+t(A))/2
  B = (B+t(B))/2
  # 3. positive definite
  if ((!check_pd(A))||(!check_pd(B))){
    stop("* measure.LERM : input is not positive definite.")
  }

  ## Main Computation
  output = tryCatch(norm(expm::logm(A)-expm::logm(B),"f"), error=function(e)e, warning=function(w)w)
  if (inherits(output, "warning")){
    message(paste("* measure.LERM :",output$message))
    return(NA)
  } else if (inherits(output, "error")){
    stop("* measure.LERM : matrix logarithm failed.")
  } else {
    return(output)
  }
}
#' @keywords internal
#' @noRd
measure.LERM.3d <- function(array3d){
  # 1. get size and set ready
  M    = dim(array3d)[3]
  logm = array(0,dim(array3d))
  # 2. a priori expm::logm computation
  for (i in 1:M){
    A = array3d[,,i]
    A = (A+t(A))/2
    output = tryCatch(expm::logm(A), error=function(e)e, warning=function(w)w)
    if (inherits(output, "warning")){
      message(paste("* measure.LERM :",output$message))
      return(NA)
    } else if (inherits(output, "error")){
      stop("* measure.LERM : matrix logarithm failed.")
    } else {
      logm[,,i] = output
    }
  }
  # 3. compute pairwise distance
  outdist = array(0,c(M,M))
  for (i in 1:(M-1)){
    logmA = logm[,,i]
    for (j in (i+1):M){
      logmB   = logm[,,j]
      normval = norm(logmA-logmB,"f")
      outdist[i,j] = normval
      outdist[j,i] = normval
    }
  }
  return(outdist)
}

# 5. KLDM -----------------------------------------------------------------
#' Symmetrized Kullback Leibler Divergence under Normal Model
#' @keywords internal
#' @noRd
measure.KLDM <- function(A,B){
  ## three conditions : size argument, symmetric, positive definite
  # 1. sqmat
  if ((!check_sqmat(A))||(!check_sqmat(B))){
    stop("* measure.KLDM : input is not square matrix/Matrix.")
  }
  # 2. symmetric
  leveltol = sqrt(.Machine$double.eps)
  if ((!isSymmetric(A, tol=leveltol))||(!isSymmetric(B, tol=leveltol))){
    stop("* measure.KLDM : input is not symmetric.")
  }
  A = (A+t(A))/2
  B = (B+t(B))/2
  # 3. positive definite
  if ((!check_pd(A))||(!check_pd(B))){
    stop("* measure.KLDM : input is not positive definite.")
  }

  ## Main Computation
  objeig = tryCatch(geigen(A, B, symmetric=TRUE, only.values = TRUE), error=function(e)e, warning=function(w)w)
  if (inherits(objeig, "error")){
    stop("* measure.KLDM : generalized eigenvalue decomposition failed.")
  } else if (inherits(objeig, "warning")){
    message(paste("* measure.KLDM :",objeig$message))
    return(NA)
  } else {
    eigvals = objeig$values
    output = sum((sqrt(eigvals)-1/sqrt(eigvals))^2)/2;
    return(output)
  }
}
#' @keywords internal
#' @noRd
measure.KLDM.3d <- function(array3d){
  # 1. get size and set ready
  M      = dim(array3d)[3]
  data3d = array(0,dim(array3d))
  # 2. a priori symmetrization
  for (i in 1:M){
    data3d[,,i] = (array3d[,,i]+t(array3d[,,i]))/2
  }
  # 3. main computation using geigen
  outdist = array(0,c(M,M))
  for (i in 1:(M-1)){
    A = data3d[,,i]
    for (j in (i+1):M){
      B = data3d[,,j]
      objeig = tryCatch(geigen(A, B, symmetric=TRUE, only.values = TRUE), error=function(e)e, warning=function(w)w)
      if (inherits(objeig, "error")){
        stop("* measure.KLDM : generalized eigenvalue decomposition failed.")
      } else if (inherits(objeig, "warning")){
        message(paste("* measure.KLDM :",objeig$message))
        return(NA)
      } else {
        eigvals = objeig$values
        outdist[i,j] = sum((sqrt(eigvals)-1/sqrt(eigvals))^2)/2;
        outdist[j,i] = sum((sqrt(eigvals)-1/sqrt(eigvals))^2)/2;
      }
    }
  }
  return(outdist)
}

# 6. JBLD -----------------------------------------------------------------
#' Jensen-Bregman Log Determinannt Divergence
#' @keywords internal
#' @noRd
measure.JBLD.3d <- function(array3d){
  # 1. get ready
  M = dim(array3d)[3]
  outdist = array(0,c(M,M))
  for (i in 1:(M-1)){
    A = array3d[,,i]
    for (j in (i+1):M){
      B = array3d[,,j]
      output = log(det((A+B)/2))-0.5*log(det(A%*%B))
      outdist[i,j] = output
      outdist[j,i] = output
    }
  }
  return(outdist)
}

# 7. Euclidean ------------------------------------------------------------
#' @keywords internal
#' @noRd
measure.Euclidean.3d <- function(array3d){
  M       = dim(array3d)[3]
  outdist = array(0, c(M,M))
  for (i in 1:(M-1)){
    A = array3d[,,i]
    for (j in (i+1):M){
      B     = array3d[,,j]
      value = norm(A-B,"f")
      outdist[i,j] = value
      outdist[j,i] = value
    }
  }
  return(outdist)
}

# 8. Cholesky -------------------------------------------------------------
#' @keywords internal
#' @noRd
measure.Choleksy.3d <- function(array3d){
  M        = dim(array3d)[3]
  outdist  = array(0, c(M,M))
  vec_chol = list()
  for (i in 1:M){
    vec_chol[[i]] = base::chol(array3d[,,i])
  }
  for (i in 1:(M-1)){
    cholA = vec_chol[[i]]
    for (j in (i+1):M){
      cholB = vec_chol[[j]]
      output = norm(cholA-cholB,"F")
      outdist[i,j] = output
      outdist[j,i] = output
    }
  }
  return(outdist)
}


# 9. RootEuclidean --------------------------------------------------------
#' @keywords internal
#' @noRd
measure.RootEuclidean.3d <- function(array3d){
  # 1. get ready for root computation
  M = dim(array3d)[3]
  root3d = array(0,dim(array3d))
  for (i in 1:M){
    root3d[,,i] = expm::sqrtm(array3d[,,i])
  }
  # 2. compute distance
  outdist = array(0,c(M,M))
  for (i in 1:(M-1)){
    A = root3d[,,i]
    for (j in (i+1):M){
      B = root3d[,,j]
      value = norm(A-B,"f")
      outdist[i,j] = value
      outdist[j,i] = value
    }
  }
  return(outdist)
}

# 10. Procrustes.SS -------------------------------------------------------
#' @keywords internal
#' @noRd
measure.Procrustes.SS.3d <- function(array3d){
  M = dim(array3d)[3]
  chol3d  = array(0,dim(array3d))
  for (i in 1:M){
    chol3d[,,i] = t(chol(array3d[,,i]))
  }
  outdist = array(0,c(M,M))
  for (i in 1:(M-1)){
    cholA = array3d[,,i]
    for (j in (i+1):M){
      cholB = array3d[,,j]
      value = pracma::procrustes(cholA,cholB)
      outdist[i,j] = value$d
      outdist[j,i] = value$d
    }
  }
  return(outdist)
}


# 11. Procrustes.Full -----------------------------------------------------
#' @keywords internal
#' @noRd
measure.Procrustes.Full.3d <- function(array3d){
  M = dim(array3d)[3]
  outdist = array(0,c(M,M))
  for (i in 1:(M-1)){
    A = array3d[,,i]
    for (j in (i+1):M){
      B = array3d[,,j]
      value = port_distcov(A,B,method="Procrustes.Full")
      outdist[i,j] = value
      outdist[j,i] = value
    }
  }
  return(outdist)
}

# 12. PowerEuclidean ------------------------------------------------------
#' @keywords internal
#' @noRd
measure.PowerEuclidean.3d <- function(A, pcoef){
  # 1.
  p = dim(A)[1]
  M = dim(A)[3]
  # 2. S_i^\alpha
  Salpha = array(0,dim(A))
  for (i in 1:M){
    eigtgt = tryCatch(eigen(A[,,i], symmetric=TRUE), error=function(e)e, warning=function(w)w)
    if (inherits(eigtgt, "error")){
      stop("* measure.PowerEuclidean : eigen decomposition failed.")
    } else if (inherits(eigtgt, "warning")){
      message(paste(" measure.PowerEuclidean :",eigtgt$warning))
      return(NA)
    } else {
      Salpha[,,i] = eigtgt$vectors %*% diag((eigtgt$values)^pcoef) %*% t(eigtgt$vectors)
    }
  }
  # 3. distance computation
  outdist = array(0,c(M,M))
  for (i in 1:(M-1)){
    A = Salpha[,,i]
    for (j in (i+1):M){
      B = Salpha[,,j]
      outvalue = norm(A-B,"f")/pcoef
      outdist[i,j] = outvalue
      outdist[j,i] = outvalue
    }
  }
}
