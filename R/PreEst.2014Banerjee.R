# Originla name : Banerjee14
#' Bayesian Estimation of a Banded Precision Matrix (Banerjee 2014)
#'
#' \code{PreEst.2014Banerjee} returns a Bayes estimator of the banded precision matrix using G-Wishart prior.
#' Stein’s loss or squared error loss function is used depending on the “loss” argument in the function.
#' The bandwidth is set at the mode of marginal posterior for the bandwidth parameter.
#'
#' @param X an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param upperK upper bound of bandwidth \eqn{k}.
#' @param delta hyperparameter for G-Wishart prior. Default value is 10. It has to be larger than 2.
#' @param logpi log of prior distribution for bandwidth \eqn{k}. Default is a function proportional to \eqn{-k^4}.
#' @param loss type of loss; either \code{"Stein"} or \code{"Squared"}.
#'
#' @return a named list containing: \describe{
#' \item{C}{a \eqn{(p\times p)} MAP estimate for precision matrix.}
#' }
#'
#' @examples
#' ## generate data from multivariate normal with Identity precision.
#' pdim = 10
#' data = matrix(rnorm(50*pdim), ncol=pdim)
#'
#' ## compare different K
#' out1 <- PreEst.2014Banerjee(data, upperK=1)
#' out2 <- PreEst.2014Banerjee(data, upperK=3)
#' out3 <- PreEst.2014Banerjee(data, upperK=5)
#'
#' ## visualize
#' opar <- par(mfrow=c(2,2), pty="s")
#' image(diag(pdim)[,pdim:1],main="Original Precision")
#' image(out1$C[,pdim:1], main="banded1::upperK=1")
#' image(out2$C[,pdim:1], main="banded1::upperK=3")
#' image(out3$C[,pdim:1], main="banded1::upperK=5")
#' par(opar)
#'
#' @references
#' \insertRef{banerjee_posterior_2014}{CovTools}
#'
#' @rdname PreEst.2014Banerjee
#' @export
PreEst.2014Banerjee <- function(X, upperK=floor(ncol(X)/2), delta=10.0, logpi=function(k){-k^4}, loss=c("Stein","Squared")){
  #-----------------------------------------------------
  ## PREPROCESSING
  #   1. typecheck : data
  fname    = "PreEst.2014Banerjee"
  checker1 = invisible_datamatrix(X, fname)
  #   2. upperK : upper bound of bandwidth
  Banerjee14.upperK = as.integer(upperK)
  if ((length(Banerjee14.upperK)!=1)||(Banerjee14.upperK<1)||(Banerjee14.upperK>ncol(X))||(abs(Banerjee14.upperK-round(Banerjee14.upperK))>sqrt(.Machine$double.eps))){
    stop("* PreEst.2014Banerjee : 'upperK' should be an integer in [1,ncol(X)].")
  }
  #   3. delta :  Default 10, has to be larger than 2
  Banerjee14.delta = as.double(delta)
  if ((length(Banerjee14.delta)!=1)||(Banerjee14.delta<=2)||(is.na(Banerjee14.delta))||(is.infinite(Banerjee14.delta))){
    stop("* PreEst.2014Banerjee : 'delta' should be a value larger than 2.")
  }
  #   4. logpi : log of prior distribution for bandwidth k.
  Banerjee14.logpi = logpi
  if (!is.function(Banerjee14.logpi)){
    stop("* PreEst.2014Banerjee : 'logpi' should be a function.")
  }
  #   5. loss : type of loss used
  Banerjee14.loss = match.arg(loss)

  #-----------------------------------------------------
  ## MAIN COMPUTATION
  output = preest.Banerjee14(Banerjee14.upperK, X, Banerjee14.delta, Banerjee14.logpi, Banerjee14.loss)

  #-----------------------------------------------------
  ## RETURN
  return(output)
}



# auxiliary functions -----------------------------------------------------
# [Banerjee14] Bayesian estimation for banded precision matrix bas --------
#######################################################################
# Bayesian estimation for banded precision matrix based on G-Wishart prior (Banerjee and Ghosal, 2014)
# K : upper bound of bandwidth
# X : n-by-p data matrix
# delta : hyperparameter for G-Wishart prior. Default value is 10. It has to be larger than 2.
# logpi : log of prior distribution for bandwidth k. Default is a function proportional to -k^4.
#
# refer : Banerjee and Ghosal, 2014. Posterior convergence rates for estimating large precision matrices using graphical models
#######################################################################
# branching case for Jay's code
#' @keywords internal
#' @noRd
preest.Banerjee14 <- function(K, X, delta, logpi, type){
  if (type=="Stein"){
    outC = preest.Banerjee14.GWishartE1(K, X, delta, logpi)
  } else if (type=="Squared"){
    outC = preest.Banerjee14.GWishartE2(K, X, delta, logpi)
  }

  output = list()
  output$C = outC
  return(output)
}
# main functions
# (1) Bayes estimator under the Stein's loss
#' @keywords internal
#' @noRd
preest.Banerjee14.GWishartE1 <- function(K, X, delta, logpi){
  X = scale(X, center = TRUE, scale =FALSE)

  logprobBG = rep(0, K)
  for(k in 1:K){
    logprobBG[k] = preest.Banerjee14.logJG(k, X, delta) + logpi(k)
  }
  khatBG = which(logprobBG == max(logprobBG)) # choice of bandwidth

  n = nrow(X); p = ncol(X)
  S = t(X)%*%X/n

  OmegaBG1 = matrix(0, p,p)
  OmegaBG1[1:(khatBG+1), 1:(khatBG+1)] = solve(diag(1/n, khatBG+1) + S[1:(khatBG+1), 1:(khatBG+1)])
  for(j in 2:(p-k)){
    OmegaBG1[j:(j+khatBG), j:(j+khatBG)] = OmegaBG1[j:(j+khatBG), j:(j+khatBG)] + solve(diag(1/n, khatBG+1) + S[j:(j+ khatBG), j:(j+khatBG)])
    OmegaBG1[j:(j+khatBG-1), j:(j+khatBG-1)] = OmegaBG1[j:(j+khatBG-1), j:(j+khatBG-1)] - solve(diag(1/n,khatBG) + S[j:(j+khatBG-1), j:(j+khatBG-1)])
  }
  OmegaBG1 = OmegaBG1*(delta+n-2)/n

  return(OmegaBG1)
}
# (2) Bayes estimator uner the Squared-error loss
#' @keywords internal
#' @noRd
preest.Banerjee14.GWishartE2 <- function(K, X, delta, logpi){
  X = scale(X, center = TRUE, scale =FALSE)

  logprobBG = rep(0, K)
  for(k in 1:K){
    logprobBG[k] = preest.Banerjee14.logJG(k, X, delta) + logpi(k)
  }
  khatBG = which(logprobBG == max(logprobBG)) # choice of bandwidth

  n = nrow(X); p = ncol(X)
  S = t(X)%*%X/n

  OmegaBG2 = matrix(0, p,p)
  OmegaBG2[1:(khatBG+1), 1:(khatBG+1)] = solve(diag(1/n,khatBG+1) + S[1:(khatBG+1), 1:(khatBG+1)])*(delta+khatBG+n)/n
  for(j in 2:(p-khatBG)){
    OmegaBG2[j:(j+khatBG), j:(j+khatBG)] = OmegaBG2[j:(j+khatBG), j:(j+khatBG)] + solve(diag(1/n,khatBG+1) + S[j:(j+khatBG), j:(j+khatBG)])*(delta+khatBG+n)/n
    OmegaBG2[j:(j+khatBG-1), j:(j+khatBG-1)] = OmegaBG2[j:(j+khatBG-1), j:(j+khatBG-1)] - solve(diag(1/n,khatBG) + S[j:(j+khatBG-1), j:(j+khatBG-1)])*(delta+khatBG+n-1)/n
  }

  return(OmegaBG2)
}
########################################################################
# auxiliary functions
########################################################################
# logarithm of J_{G^k}(delta, n, I_p, nS) in Banerjee and Ghosal (2014) (p2123: eq (5.7))
#' @keywords internal
#' @noRd
preest.Banerjee14.logJG <- function(k, X, delta){
  n = nrow(X); p = ncol(X)
  S = t(X)%*%X/n

  res = 0
  res = res + sum(lgamma(0.5*(delta+n+(0:k))) - lgamma(0.5*(delta+(0:k))))
  res = res + (p-k-1)*(lgamma(0.5*(delta+n+k)) - lgamma(0.5*(delta+k)))
  SC1 = S[1:(1+k), 1:(1+k)]
  res = res - 0.5*(delta+n+k)*determinant(diag(1,nrow=k+1)+n*SC1, logarithm = T)$modulus
  for(j in 2:(p-k)){
    SSj = S[j:(j+k-1), j:(j+k-1)]
    SCj = S[j:(j+k), j:(j+k)]
    res = res + 0.5*(delta+n+k-1)*determinant(diag(1,nrow=k)+n*SSj, logarithm = T)$modulus - 0.5*(delta+n+k)*determinant(diag(1,nrow=k+1)+n*SCj, logarithm = T)$modulus
  }

  return(res)
}

