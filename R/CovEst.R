#' Estimating Covariance Matrix
#'
#' It compiles several popular methods of inferring covariance structure from observed data.
#' Several principles or assumed structures and corresponding estimation techniques are included,
#' such as sparse covariance estimation (\code{Bickel08}, \code{Cai11}), soft thresholding (\code{Donoho95}),
#' near positive definiteness (\code{Qi06}), and so on.
#'
#' @param X an \code{(n-by-p)} matrix where each row is an observation from the first dataset.
#' @param method a name of estimation method.
#' @param param a parameter value >= 0.
#'
#'
#' @return a named list containing: \describe{
#' \item{S}{a \code{(p-by-p)} covariance matrix estimate.}
#' }
#'
#' @examples
#' ## generate data from multivariate normal with Identity covariance.
#' data = mvtnorm::rmvnorm(100, sigma=diag(10))
#'
#' ## run estimation
#' out3 = CovEst(data, method="Donoho95")
#' out4 = CovEst(data, method="Fan13")
#' out5 = CovEst(data, method="Qi06")
#'
#' ## Visualize
#' par(mfrow=c(2,3))
#' image(pracma::flipud(diag(10)),main="Original Covariance")
#' image(pracma::flipud(out3$S), main="Donoho95")
#' image(pracma::flipud(out4$S), main="Fan13")
#' image(pracma::flipud(out5$S), main="Qi06")
#'
#' @references [Donoho95] Donoho, D. et al. (1995) \emph{Wavelet Shrinkage: Asymptopia?} Journal of the Royal Statistical Society Series B, Vol.57(2):301-369.
#' @references [Fan13] Fan. J. et al. (2013) \emph{Large covariance estimation by thresholding principal orthogonal complements.}
#' Journal of the Royal Statistical Society Series B, Vol.75(4):603-680.
#' @references [Qi06] Qi, H. and Sun, D. (2006) \emph{A Quadratically Convergent Newton Method for Computing the Nearest Correlation Matrix.} SIAM J.Matrix Anal.& Appl., Vol.28(2):360-385.
#' @export
CovEst <- function(X, method=c("Bickel08","Cai11","Donoho95","Fan13","Qi06"),param=1.0){
  ## PREPROCESSING ON INPUTS AND PARAMETERS
  # 1. valid data matrix
  if (!check_datamatrix(X)){stop("* CovEst : an input data matrix X is invalid.")}
  if ((nrow(X)==1)||(ncol(X)==1)){stop("* CovEst : invalid input matrix X.")}
  # 2. method : THIS SHOULD BE UPDATED EVERYTIME A METHOD IS ADDED
  if (missing(method)){method="Bickel08"} else {
    method = match.arg(method)
  }

  ## Parameter Checking
  if ((length(param)!=1)||(param<0)||(is.infinite(param))||(!is.numeric(param))||(is.na(param))){
    stop("* CovEst : 'param' value should be >=0.")
  }

  ## Main Computation
  output = switch(method,
                  Donoho95 = covest.Donoho95.once(X,param),  # threshold of soft thresholding
                  Fan13    = covest.Fan13.once(X,param),
                  Qi06     = covest.Qi06.once(X)             # anyway, no CV originally.
  )
  ## RETURN OUTPUT
  return(output)
}
