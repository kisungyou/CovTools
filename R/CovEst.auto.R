#' Estimating Covariance Matrix with automatic tuning
#'
#' In this variant of \code{\link{CovEst}}, we implemented automatic parameter tuning scheme
#' applying 2-fold cross validation repeatedly and choosing the minimal one with the least discrepancy.
#'
#' @param X an \code{(n-by-p)} matrix where each row is an observation from the first dataset.
#' @param method a name of estimation method.
#' @param opt a list of options containing following fields: \describe{
#' \item{nCV}{the number for repetitions for 2-fold random cross validation.}
#' \item{nsearch}{the number of trials on range of regularization parameters.}
#' }
#' @param parallel a logical; \code{TRUE} to use half of available cores, \code{FALSE} to do every computation sequentially.
#'
#'
#' @return a named list containing: \describe{
#' \item{S}{a \code{(p-by-p)} covariance matrix estimate.}
#' \item{CV}{a dataframe containing vector of tested threshold values(\code{thr}) and corresponding cross validation scores(\code{CVscore}).}
#' }
#'
#' @examples
#' \dontrun{
#' ## generate data from multivariate normal with Identity covariance.
#' data = mvtnorm::rmvnorm(100, sigma=diag(5))
#'
#' ## run automatic estimation
#' sopt = list(nCV=2,nsearch=3) # common option
#' out1 = CovEst.auto(data, method="Bickel08", opt=sopt)
#' out2 = CovEst.auto(data, method="Cai11",    opt=sopt)
#'
#' ## Visualize
#' par(mfrow=c(1,3))
#' image(pracma::flipud(diag(5)),main="Original Covariance")
#' image(pracma::flipud(out1$S), main="Bickel08")
#' image(pracma::flipud(out2$S), main="Cai11")

#' }
#'
#' @references [Cai11] Cai, T. and Liu, W. (2011) \emph{Adaptive Thresholding for Sparse Covariance Matrix Estimation.} Journal of the American Statistical Association, Vol.106:672-684.
#' @references [Donoho95] Donoho, D. et al. (1995) \emph{Wavelet Shrinkage: Asymptopia?} Journal of the Royal Statistical Society Series B, Vol.57(2):301-369.
#' @references [Fan13] Fan. J. et al. (2013) \emph{Large covariance estimation by thresholding principal orthogonal complements.}
#' Journal of the Royal Statistical Society Series B, Vol.75(4):603-680.
#' @references [Qi06] Qi, H. and Sun, D. (2006) \emph{A Quadratically Convergent Newton Method for Computing the Nearest Correlation Matrix.} SIAM J.Matrix Anal.& Appl., Vol.28(2):360-385.
#'
#'
#' @seealso \code{\link{CovEst}}
#' @export
CovEst.auto <- function(X, method=c("Bickel08","Cai11","Donoho95","Fan13","Qi06"),opt=list(nCV=10,nsearch=10),parallel=FALSE){
  ## PREPROCESSING ON INPUTS AND PARAMETERS
  # 1. valid data matrix
  if (!check_datamatrix(X)){stop("* CovEst.auto : an input data matrix X is invalid.")}
  if ((nrow(X)==1)||(ncol(X)==1)){stop("* CovEst.auto : invalid input matrix X.")}
  # 2. method : THIS SHOULD BE UPDATED EVERYTIME A METHOD IS ADDED
  if (missing(method)){method="Bickel08"} else {
    method = match.arg(method)
  }

  ## Parameter Checking from Option
  if (!missing(opt)){
    # 1. nCV : cross validation numerics
    if ("nCV" %in% names(opt)){
      nCV = opt$nCV
      if ((length(nCV)!=1)||(!is.numeric(nCV))||(nCV<1)||(nCV>nrow(X))||is.na(nCV)||(abs(nCV-round(nCV))>sqrt(.Machine$double.eps))){
        stop("* CovEst.auto : opt$nCV parameter is not valid.")
      }
      nCV = round(nCV)
    } else {nCV = round(10) }
    # 2. nsearch : grid generation
    if ("nsearch" %in% names(opt)){
      nsearch = opt$nsearch
      if ((length(nsearch)!=1)||(!is.numeric(nsearch))||(nsearch<1)||(nsearch>nrow(X))||is.na(nsearch)||(abs(nsearch-round(nsearch))>sqrt(.Machine$double.eps))){
        stop("* CovEst.auto : opt$nsearch parameter is not valid.")
      }
      nsearch = round(nsearch)
    } else {
      nsearch = round(10)
    }
  } else { # missing cases
    nCV = round(10)
    nsearch = round(10)
  }





  ## Parallel Setting
  if (!parallel){nCore = 1  }
  else {nCore = max(round(detectCores()/2),1)}

  ## Main Computation
  output = switch(method,
                  Donoho95 = covest.Donoho95(X,nCV,nCore,nsearch),
                  Fan13    = covest.Fan13(X,nCV,nCore,nsearch),
                  Qi06     = covest.Qi06(X)
  )
  ## RETURN OUTPUT
  return(output)
}
