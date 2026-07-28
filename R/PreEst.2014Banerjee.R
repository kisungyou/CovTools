# Originla name : Banerjee14
#' Bayesian Estimation of a Banded Precision Matrix (Banerjee 2014)
#'
#' \code{PreEst.2014Banerjee} returns a Bayes estimator of the banded precision matrix using G-Wishart prior.
#' Stein’s loss or squared error loss function is used depending on the “loss” argument in the function.
#' The bandwidth is set at the mode of marginal posterior for the bandwidth parameter.
#'
#' @param X an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param upperK integer upper bound of bandwidth \eqn{k}, at most \eqn{p-1}.
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
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
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
  if (!is.numeric(X)){
    stop("* PreEst.2014Banerjee : 'X' should be a numeric matrix.")
  }
  X = as.matrix(X)
  n = nrow(X)
  p = ncol(X)
  if ((n < 2)||(p < 2)){
    stop("* PreEst.2014Banerjee : 'X' should have at least 2 rows and 2 columns.")
  }
  #   2. upperK : upper bound of bandwidth
  if ((!is.numeric(upperK))||(length(upperK)!=1)||(!is.finite(upperK))||
      (upperK != floor(upperK))||(upperK<1)||(upperK>(p-1))){
    stop("* PreEst.2014Banerjee : 'upperK' should be an integer in [1,ncol(X)-1].")
  }
  Banerjee14.upperK = as.integer(upperK)
  #   3. delta :  Default 10, has to be larger than 2
  if ((!is.numeric(delta))||(length(delta)!=1)||(!is.finite(delta))||(delta<=2)){
    stop("* PreEst.2014Banerjee : 'delta' should be a value larger than 2.")
  }
  Banerjee14.delta = as.double(delta)
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

  khatBG = preest.Banerjee14.selectK(K, X, delta, logpi)

  n = nrow(X); p = ncol(X)
  S = t(X)%*%X/n

  OmegaBG1 = matrix(0, p,p)
  OmegaBG1[1:(khatBG+1), 1:(khatBG+1)] = solve(diag(1/n, khatBG+1) + S[1:(khatBG+1), 1:(khatBG+1)])
  if ((p-khatBG) >= 2){
    for(j in seq.int(2, p-khatBG)){
      OmegaBG1[j:(j+khatBG), j:(j+khatBG)] = OmegaBG1[j:(j+khatBG), j:(j+khatBG)] + solve(diag(1/n, khatBG+1) + S[j:(j+ khatBG), j:(j+khatBG)])
      OmegaBG1[j:(j+khatBG-1), j:(j+khatBG-1)] = OmegaBG1[j:(j+khatBG-1), j:(j+khatBG-1)] - solve(diag(1/n,khatBG) + S[j:(j+khatBG-1), j:(j+khatBG-1)])
    }
  }
  OmegaBG1 = OmegaBG1*(delta+n-2)/n

  return(preest.Banerjee14.ensurePD(OmegaBG1))
}
# (2) Bayes estimator uner the Squared-error loss
#' @keywords internal
#' @noRd
preest.Banerjee14.GWishartE2 <- function(K, X, delta, logpi){
  X = scale(X, center = TRUE, scale =FALSE)

  khatBG = preest.Banerjee14.selectK(K, X, delta, logpi)

  n = nrow(X); p = ncol(X)
  S = t(X)%*%X/n

  OmegaBG2 = matrix(0, p,p)
  OmegaBG2[1:(khatBG+1), 1:(khatBG+1)] = solve(diag(1/n,khatBG+1) + S[1:(khatBG+1), 1:(khatBG+1)])*(delta+khatBG+n)/n
  if ((p-khatBG) >= 2){
    for(j in seq.int(2, p-khatBG)){
      OmegaBG2[j:(j+khatBG), j:(j+khatBG)] = OmegaBG2[j:(j+khatBG), j:(j+khatBG)] + solve(diag(1/n,khatBG+1) + S[j:(j+khatBG), j:(j+khatBG)])*(delta+khatBG+n)/n
      OmegaBG2[j:(j+khatBG-1), j:(j+khatBG-1)] = OmegaBG2[j:(j+khatBG-1), j:(j+khatBG-1)] - solve(diag(1/n,khatBG) + S[j:(j+khatBG-1), j:(j+khatBG-1)])*(delta+khatBG+n-1)/n
    }
  }

  return(preest.Banerjee14.ensurePD(OmegaBG2))
}

# Select a single posterior mode. which.max() intentionally resolves ties by
# choosing the smallest bandwidth rather than returning a vector index.
#' @keywords internal
#' @noRd
preest.Banerjee14.selectK <- function(K, X, delta, logpi){
  logprobBG = vapply(seq_len(K), function(k){
    logprior = logpi(k)
    if ((!is.numeric(logprior))||(length(logprior)!=1)||(!is.finite(logprior))){
      stop("* PreEst.2014Banerjee : 'logpi' should return one finite numeric value.")
    }
    score = as.numeric(preest.Banerjee14.logJG(k, X, delta)) + logprior
    if (!is.finite(score)){
      stop("* PreEst.2014Banerjee : posterior bandwidth score is not finite.")
    }
    score
  }, numeric(1))
  best = max(logprobBG)
  tolerance = 1000*.Machine$double.eps*max(1, max(abs(logprobBG)))
  return(as.integer(which(logprobBG >= (best-tolerance))[1]))
}

#' @keywords internal
#' @noRd
preest.Banerjee14.ensurePD <- function(Omega){
  Omega = (Omega+t(Omega))/2
  if (any(!is.finite(Omega))){
    stop("* PreEst.2014Banerjee : precision estimate contains non-finite values.")
  }
  okay = tryCatch({
    chol(Omega)
    TRUE
  }, error=function(e) FALSE)
  if (!okay){
    stop("* PreEst.2014Banerjee : precision estimate is not positive definite.")
  }
  return(Omega)
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
  if ((p-k) >= 2){
    for(j in seq.int(2, p-k)){
      SSj = S[j:(j+k-1), j:(j+k-1)]
      SCj = S[j:(j+k), j:(j+k)]
      res = res + 0.5*(delta+n+k-1)*determinant(diag(1,nrow=k)+n*SSj, logarithm = TRUE)$modulus - 0.5*(delta+n+k)*determinant(diag(1,nrow=k+1)+n*SCj, logarithm = TRUE)$modulus
    }
  }

  return(res)
}
