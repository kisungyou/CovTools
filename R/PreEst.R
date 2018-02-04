#' Estimating Precision Matrix
#'
#' This code compiles several estimation methods for \emph{precision} matrix, which is
#' an inverse of covariance matrix, including penalized likelihood method with L1 penalty(\code{Yuan07}, \code{Banerjee06})
#' or Bayesian approaches incorporating banded structure assumptions(\code{Banerjee14},\code{Lee17}).
#'
#' @param X an \code{(n-by-p)} matrix where each row is an observation from the first dataset.
#' @param method a name of estimation method.
#' @param opt a list containing following parameters, \describe{
#' \item{Banerjee06.confidence}{level of confidence in \eqn{(0,1)}.}
#' \item{Banerjee14.upperK}{upper bound for bandwidth.}
#' \item{Banerjee14.delta}{a number larger than 2.}
#' \item{Banerjee14.logpi}{log of prior distribution for bandwidth \code{k}.}
#' \item{Banerjee14.loss}{loss type; either \code{"Stein"} or \code{"Squared"} type.}
#' \item{Lee17.upperK}{upper bound for bandwidth.}
#' \item{Lee17.logpi}{log of prior distribution for bandwidth \code{k}.}
#' \item{Yuan07.lambda}{a regularization parameter.}
#' }
#'
#' @return a \code{(p-by-p)} estimated precision matrix.
#'
#' @examples
#' ## generate data from multivariate normal with Identity precision.
#' data = mvtnorm::rmvnorm(100, sigma=diag(10))
#'
#' ## run estimation
#' out1 = PreEst(data, method="Banerjee06")
#' out2 = PreEst(data, method="Banerjee14")
#' out3 = PreEst(data, method="Lee17")
#' out4 = PreEst(data, method="Yuan07")
#'
#' ## Visualize
#' par(mfrow=c(2,2))
#' image(pracma::flipud(out1), main="Banerjee06")
#' image(pracma::flipud(out2), main="Banerjee14")
#' image(pracma::flipud(out3), main="Lee17")
#' image(pracma::flipud(out4), main="Yuan07")
#'
#' @references [Banerjee06] Banerjee et al (2006) \emph{Convex optimization techniques for fitting sparse Gaussian graphical models.} ICML'06:89-96.
#' @references [Banerjee14] Banerjee, S. and Ghosal, S. (2014) \emph{Posterior convergence rates for estimating large precision matrices using graphical models.} Electronic Journal of Statistics, Vol.8:2111-2137.
#' @references [Lee17] Lee, K. and Lee, J. (2017) \emph{Estimating Large Precision Matrices via Modified Cholesky Decomposition.} arXiv:1707.01143.
#' @references [Yuan07] Yuan, M. and Lin, Y. (2007) \emph{Model selection and estimation in the Gaussian graphical model.} Biometrika, Vol.94(1):19-35.
#'
#' @export
PreEst<- function(X, method=c("Banerjee06","Banerjee14","Lee17","Yuan07"),
                  opt=list(Banerjee06.confidence=0.95, # this part should be 'lambda' for nonauto.
                           Banerjee14.upperK = floor(ncol(X)/2),
                           Banerjee14.delta = 10.0,
                           Banerjee14.logpi=function(k){-k^4},
                           Banerjee14.loss=c("Stein","Squared"),
                           Lee17.upperK=floor(ncol(X)/2),
                           Lee17.logpi=function(k){-k^4},
                           Yuan07.lambda=1.0)){

  ## PREPROCESSING ON INPUTS AND PARAMETERS
  # 1. valid data matrix
  if (!check_datamatrix(X)){stop("* PreEst : an input data matrix X is invalid.")}
  if ((nrow(X)==1)||(ncol(X)==1)){stop("* PreEst : invalid input matrix X.")}
  # 2. method : THIS SHOULD BE UPDATED EVERYTIME A METHOD IS ADDED
  if (missing(method)){method="Yuan07"} else {
    method = match.arg(method)
  }

  ## Parameter Checking from Option
  if (!missing(opt)){
    # 1. Yuan07.lambda :
    if ("Yuan07.lambda" %in% names(opt)){
      Yuan07.lambda = opt$Yuan07.lambda
      if ((length(Yuan07.lambda)!=1)||(is.na(Yuan07.lambda))||(is.infinite(Yuan07.lambda))||(Yuan07.lambda<0)){
        stop("* PreEst : 'Yuan07.lambda' should be a positive real number.")
      }
      Yuan07.lambda = as.double(Yuan07.lambda)
    } else {
      Yuan07.lambda = 1.0
    }
    # 2. Banerjee06.confidence : for automatic selection (0,1)
    if ("Banerjee06.confidence" %in% names(opt)){
      Banerjee06.confidence = opt$Banerjee06.confidence
      if ((length(Banerjee06.confidence)!=1)||(Banerjee06.confidence<=0)||(Banerjee06.confidence>=1)||(is.na(Banerjee06.confidence))){
        stop("* PreEst : 'Banerjee06.confidence' should be a real number in (0,1).")
      }
    } else {
      Banerjee06.confidence = 0.95
    }
    # 3. Lee17.upperK : upper bound of bandwidth
    if ("Lee17.upperK" %in% names(opt)){
      Lee17.upperK = opt$Lee17.upperK
      if ((length(Lee17.upperK)!=1)||(Lee17.upperK<1)||(Lee17.upperK>ncol(X))||(abs(Lee17.upperK-round(Lee17.upperK))>sqrt(.Machine$double.eps))){
        stop("* PreEst : 'Lee17.upperK' should be an integer in [1,ncol(X)].")
      }
      Lee17.upperK = round(Lee17.upperK)
    } else {
      Lee17.upperK=floor(ncol(X)/2)
    }
    # 4. Lee17.logpi : log of prior distribution for bandwidth k.
    if ("Lee17.logpi" %in% names(opt)){
      Lee17.logpi = opt$Lee17.logpi
      if (!is.function(Lee17.logpi)){
        stop("* PreEst : 'Lee17.logpi' should be a function.")
      }
    } else {
      Lee17.logpi = function(k){-k^4}
    }
    # 5. Banerjee14.upperK
    if ("Banerjee14.upperK" %in% names(opt)){
      Banerjee14.upperK = opt$Banerjee14.upperK
      if ((length(Banerjee14.upperK)!=1)||(Banerjee14.upperK<1)||(Banerjee14.upperK>ncol(X))||(abs(Banerjee14.upperK-round(Banerjee14.upperK))>sqrt(.Machine$double.eps))){
        stop("* PreEst : 'Banerjee14.upperK' should be an integer in [1,ncol(X)].")
      }
      Banerjee14.upperK = round(Banerjee14.upperK)
    } else {
      Banerjee14.upperK=floor(ncol(X)/2)
    }
    # 6. Banerjee14.delta : Default 10, has to be larger than 2
    if ("Banerjee14.delta" %in% names(opt)){
      Banerjee14.delta = opt$Banerjee14.delta
      if ((length(Banerjee14.delta)!=1)||(Banerjee14.delta<=2)||(is.na(Banerjee14.delta))||(is.infinite(Banerjee14.delta))){
        stop("* PreEst : 'Banerjee14.delta' should be a value larger than 2.")
      }
    } else {
      Banerjee14.delta = 10.0
    }
    # 7. Banerjee14.logpi : same as Lee17.logpi
    if ("Banerjee14.logpi" %in% names(opt)){
      Banerjee14.logpi = opt$Banerjee14.logpi
      if (!is.function(Banerjee14.logpi)){
        stop("* PreEst : 'Banerjee14.logpi' should be a function.")
      }
    } else {
      Banerjee14.logpi = function(k){-k^4}
    }
    # 8. Banerjee14.loss : type of loss used
    if ("Banerjee14.loss" %in% names(opt)){
      Banerjee14.loss = opt$Banerjee14.loss
      if (!(Banerjee14.loss %in% c("Stein","Squared"))){
        stop("* PreEst : 'Banerjee14.loss' should be either 'Stein' or 'Squared'")
      }
    } else {
      Banerjee14.loss = "Stein"
    }
  } else { # if missing
    Yuan07.lambda = 1.0
    Banerjee06.confidence = 0.95
    Lee17.upperK=floor(ncol(X)/2)
    Lee17.logpi=function(k){-k^4}
    Banerjee14.upperK=floor(ncol(X)/2)
    Banerjee14.delta=10.0
    Banerjee14.logpi = function(k){-k^4}
    Banerjee14.loss  = "Stein"
  }

  ## Main Computation
  output = switch(method,
                  Yuan07     = preest.Yuan07.once.matrix(X, Yuan07.lambda),
                  Lee17      = preest.Lee17(Lee17.upperK, X, logpi=Lee17.logpi),
                  Banerjee06 = preest.Banerjee06(X, Banerjee06.confidence),
                  Banerjee14 = preest.Banerjee14(Banerjee14.upperK, X, Banerjee14.delta, Banerjee14.logpi, Banerjee14.loss)
  )
  ## RETURN OUTPUT
  return(output$C)
}
