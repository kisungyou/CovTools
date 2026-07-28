# Original Article : An.2014 + Bickel_Levina.2008
#' Banded Precision Matrix Estimation via Bandwidth Test
#'
#' \code{PreEst.2014An} returns an estimator of the banded precision matrix using the modified Cholesky decomposition.
#' It uses the estimator defined in Bickel and Levina (2008). The bandwidth is determined by the bandwidth test
#' suggested by An, Guo and Liu (2014).
#'
#' @param X an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param upperK integer upper bound of bandwidth \eqn{k}, at most
#' \eqn{\min(p-1,n-5)}.
#' @param algorithm bandwidth test algorithm to be used.
#' @param alpha significance level for the bandwidth test.
#'
#' @return a named list containing: \describe{
#' \item{C}{a \eqn{(p\times p)} estimated banded precision matrix.}
#' \item{optk}{an estimated optimal bandwidth acquired from the test procedure.}
#' }
#'
#' @examples
#' \dontrun{
#' ## parameter setting
#' p = 200; n = 100
#' k0 = 5; A0min=0.1; A0max=0.2; D0min=2; D0max=5
#'
#' set.seed(123)
#' A0 = matrix(0, p,p)
#' for(i in 2:p){
#'   term1 = runif(n=min(k0,i-1),min=A0min, max=A0max)
#'   term2 = sample(c(1,-1),size=min(k0,i-1),replace=TRUE)
#'   vals  = term1*term2
#'   vals  = vals[ order(abs(vals)) ]
#'   A0[i, max(1, i-k0):(i-1)] = vals
#' }
#'
#' D0 = diag(runif(n = p, min = D0min, max = D0max))
#' Omega0 = t(diag(p) - A0)%*%diag(1/diag(D0))%*%(diag(p) - A0)
#'
#' ## data generation (based on AR representation)
#' ## it is same with generating n random samples from N_p(0, Omega0^{-1})
#' X = matrix(0, nrow=n, ncol=p)
#' X[,1] = rnorm(n, sd = sqrt(D0[1,1]))
#' for(j in 2:p){
#'   mean.vec.j = X[, 1:(j-1)]%*%as.matrix(A0[j, 1:(j-1)])
#'   X[,j] = rnorm(n, mean = mean.vec.j, sd = sqrt(D0[j,j]))
#' }
#'
#' ## banded estimation using two different schemes
#' Omega1 <- PreEst.2014An(X, upperK=20, algorithm="Bonferroni")
#' Omega2 <- PreEst.2014An(X, upperK=20, algorithm="Holm")
#'
#' ## visualize true and estimated precision matrices
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' image(Omega0[,p:1],   main="Original Precision")
#' image(Omega1$C[,p:1], main="banded3::Bonferroni")
#' image(Omega2$C[,p:1], main="banded3::Holm")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{an_hypothesis_2014}{CovTools}
#'
#' \insertRef{bickel_regularized_2008}{CovTools}
#'
#' @rdname PreEst.2014An
#' @export
PreEst.2014An <- function(X, upperK=min(floor(ncol(X)/2), ncol(X)-1, nrow(X)-5),
                          algorithm=c("Bonferroni","Holm"), alpha=0.01){
  #-----------------------------------------------------
  ## PREPROCESSING
  #   1. typecheck : data
  fname    = "PreEst.2014An"
  checker1 = invisible_datamatrix(X, fname)
  if (!is.numeric(X)){
    stop("* PreEst.2014An : 'X' should be a numeric matrix.")
  }
  X = as.matrix(X)
  n = nrow(X)
  p = ncol(X)
  if ((n < 6)||(p < 2)){
    stop("* PreEst.2014An : 'X' should have at least 6 rows and 2 columns.")
  }
  #   2. upperK
  maxK = min(p-1, n-5)
  if ((!is.numeric(upperK))||(length(upperK)!=1)||(!is.finite(upperK))||
      (upperK != floor(upperK))||(upperK<1)||(upperK>maxK)){
    stop(paste0("* PreEst.2014An : 'upperK' should be an integer in [1,",
                maxK, "]."))
  }
  upperK = as.integer(upperK)
  #   3. algorithm
  if (missing(algorithm)){
    algorithm = "Bonferroni"
  } else {
    algorithm = match.arg(algorithm)
  }
  #   4. alpha
  if ((!is.numeric(alpha))||(length(alpha)!=1)||(!is.finite(alpha))||
      (alpha<=0)||(alpha>=1)){
    stop("* PreEst.2014An : 'alpha' should be a real number in (0,1).")
  }
  alpha = as.double(alpha)

  # Precision estimation is invariant to the unknown mean.
  X = scale(X, center=TRUE, scale=FALSE)
  if (any(colSums(X^2) <= .Machine$double.eps)){
    stop("* PreEst.2014An : every column of 'X' should have positive variance.")
  }

  #-----------------------------------------------------
  ##  MAIN COMPUTATION
  # An's method
  K = upperK
  pk.vec = rep(0, K)
  for (k in 1:K){
    pk.vec[k] = p.bandAN(k, X)
  }
  if (any(!is.finite(pk.vec))){
    stop("* PreEst.2014An : bandwidth test produced a non-finite p-value.")
  }

  # bandwidth selection procedure
  if(algorithm=="Bonferroni"){ # algorithm 1
    rejected = which(pk.vec <= alpha/K)
    khat = if (length(rejected)==0) 1L else max(rejected)
  }else{ # algorithm 2
    ordering = order(pk.vec)
    ordered.p = pk.vec[ordering]
    alpha.vec = alpha/(K-seq_len(K)+1)
    first.not.rejected = which(ordered.p > alpha.vec)[1]
    nrejected = if (is.na(first.not.rejected)) K else first.not.rejected-1L
    khat = if (nrejected==0) 1L else max(ordering[seq_len(nrejected)])
  }
  khat = as.integer(khat)

  # Estimating the precision matrix based on khat
  Dhat = diag(p)
  Ahat = matrix(0, p,p)
  for(j in 1:p){
    dvalue = as.numeric(dhat(j, khat, X))
    if ((length(dvalue)!=1)||(!is.finite(dvalue))||(dvalue<=0)){
      stop("* PreEst.2014An : residual variance estimate is not positive.")
    }
    Dhat[j,j] = dvalue # ordinary LSE
  }
  for(j in 2:p){
    Zj = Z(j, khat, X)
    Xj = X[, j]
    Covhat = t(Zj)%*%matrix(Xj)/nrow(X)   # k X 1 matrix
    ahat = An14.solve(VhatZ(j, khat, X), Covhat) # k X 1 matrix
    Ahat[j, max(1, j-khat):(j-1)] = c(ahat)
  }
  Omegahat = t(diag(p) - Ahat)%*%diag(1/diag(Dhat))%*%(diag(p) - Ahat) # estimated precision matrix
  Omegahat = (Omegahat+t(Omegahat))/2

  #-----------------------------------------------------
  ## RETURN OUTPUT
  output = list()
  output$C = Omegahat
  output$optk = khat
  return(output)
}





#######################################################################
# Auxiliary functions for band test in An et al. (2014)
#######################################################################
#' @keywords internal
#' @noRd
p.bandAN <- function(k, X){
  res = 2 * (1 - pnorm( abs(LfAN(k, X)) ))
  return(res)
}
#' @keywords internal
#' @noRd
# Z_j(k) : n X k matrix
Z <- function(j, k, X){
  ind = c(max(1, j-k):(j-1))
  res = X[, ind]
  if(j == 2 || k == 1) res = matrix(res)
  return(res)
}
#' @keywords internal
#' @noRd
# \hat{Var}(Z_j(k)) : k X k matrix
VhatZ <- function(j, k, X){
  Zj = Z(j, k, X)
  res = t(Zj)%*%Zj/nrow(X)
  return(res)
}
#' @keywords internal
#' @noRd
# \hat{d}_{jk}
dhat <- function(j, k, X){
  n = nrow(X)
  Xj = X[, j]
  if(j == 1){
    res = sum(Xj^2)/n   # scalar
  }else{
    Zj = Z(j, k, X)
    VhatX = sum(Xj^2)/n   # scalar
    Covhat = t(Zj)%*%matrix(Xj)/n   # k X 1 matrix
    # res = VhatX - t(Covhat)%*%solve( VhatZ(j, k, X), Covhat)
    res = VhatX - t(Covhat)%*%An14.solve(VhatZ(j, k, X), Covhat)
    if(res <= 0){
      ahat = An14.solve(VhatZ(j, k, X), Covhat)
      res = t(as.matrix(c(-ahat, 1))) %*% VhatZ(j+1,k+1,X) %*% as.matrix(c(-ahat, 1))
    }
  }
  return(res)
}
#' @keywords internal
#' @noRd
# \hat{d}_{jk} in An et al. (2014)
dhatAN <- function(j, k, X){
  n = nrow(X)
  Xj = X[, j]
  if(j == 1){
    res = sum(Xj^2)   # scalar
  }else{
    Zj = Z(j, k, X)
    nVhatX = sum(Xj^2)   # scalar
    nCovhat = t(Zj)%*%matrix(Xj)   # k X 1 matrix
    res = nVhatX - t(nCovhat)%*%An14.solve(n*VhatZ(j, k, X), nCovhat)
  }
  return(res/(n - j + max(1,j-k)))
}
#' @keywords internal
#' @noRd
thatAN <- function(j, k, X){ # scalar
  n = nrow(X)
  Xj = X[, j]
  Zj = Z(j, k, X)

  Covhat = t(Zj)%*%matrix(Xj)/n   # k X 1 matrix
  res = An14.solve(VhatZ(j, k, X), Covhat)
  return(res[1])
}
#' @keywords internal
#' @noRd
inv.deltahatAN <- function(j, k, X){
  n = nrow(X)

  Xjk = matrix(X[, j-k])
  if(k == 1){
    res = sum(Xjk^2)
  }else{
    Zjk1 = Z(j, k-1, X)
    res = (sum(Xjk^2) -
      t(Xjk)%*%Zjk1%*%An14.solve(t(Zjk1)%*%Zjk1, t(Zjk1)%*%Xjk))
  }
  return(res)
}
#' @keywords internal
#' @noRd
LfAN <- function(k, X){
  n = nrow(X)
  p = ncol(X)

  Lf = 0
  for(j in (k+1):p){
    Lf = Lf + thatAN(j, k, X)^2 / dhatAN(j, k, X) * inv.deltahatAN(j, k, X) - (n-k)/(n-k-2)
  }
  Lf = Lf * (n-k-2)*sqrt(n-k-4)/(n-k)/sqrt(2*(n-k-1))/sqrt(p-k)
  return(as.numeric(Lf))
}

#' @keywords internal
#' @noRd
An14.solve <- function(A, B=NULL){
  output = tryCatch({
    if (is.null(B)){
      solve(A)
    } else {
      solve(A, B)
    }
  }, error=function(e){
    stop("* PreEst.2014An : a regression design is singular.", call.=FALSE)
  })
  if (any(!is.finite(output))){
    stop("* PreEst.2014An : a regression solve produced non-finite values.")
  }
  return(output)
}
