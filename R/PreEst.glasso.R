#' Precision Matrix Estimation via Graphical Lasso
#'
#' Given a sample covariance matrix \eqn{S}, graphical lasso aims at estimating sparse precision matrix \eqn{X} - inverse
#' of covariance. It solves a following optimization problem,
#' \deqn{\textrm{max}_X
#' \log\textrm{det}X - <S,X> - \lambda \|X \|_1 \textrm{ such that } X \succ 0}
#' where \eqn{\lambda} a regularization parameter, \eqn{<S,X>=tr(S^T X)} , \eqn{\|X\|_1 = \sum X_{ij}} and \eqn{X\succ 0} indicates positive definiteness. We provide three
#' modes of computations, \code{'fixed'},\code{'confidence'}, or \code{'BIC'} with respect to \eqn{\lambda}. Please see the section below for more details.
#'
#' @section regularization parameters:
#' We currently provide three options for solving the problem, \code{'fixed'},\code{'confidence'}, or \code{'BIC'} with respect to \eqn{\lambda}.
#' When the method type is \code{'fixed'}, the parameter should be a single numeric value as a user-defined \eqn{\lambda} value. Likewise,
#' method type of \code{'confidence'} requires a singule numeric value in \eqn{(0,1)}, where the value is set heuristically
#' according to
#' \deqn{
#' \rho = \frac{t_{n-2}(\gamma) \max S_{ii}S_{jj}}{\sqrt{n-2+ t_{n-2}^2(\gamma)}}
#' }
#' for a given confidence level \eqn{\gamma \in (0,1)} as proposed by Banerjee et al. (2006).
#' Finally, \code{'BIC'} type requires a vector of \eqn{\lambda} values and opts for a lambda value with the lowest BIC values
#' as proposed by Yuan and Lin (2007).
#'
#'
#' @param X an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param method a list containing following parameters, \describe{
#' \item{type}{one of \code{'fixed'},\code{'confidence'}, or \code{'BIC'}.}
#' \item{param}{either a numeric value or vector of values.}
#' }
#' @param parallel a logical; \code{TRUE} for using half the cores available, \code{FALSE} otherwise.
#'
#' @return a named list containing: \describe{
#' \item{C}{a \eqn{(p\times p)} estimated precision matrix.}
#' \item{BIC}{a dataframe containing \eqn{\lambda} values and corresponding BIC scores with \code{type='BIC'} method.}
#' }
#'
#' @examples
#' ## generate data from multivariate normal with Identity precision.
#' data = mvtnorm::rmvnorm(100, sigma=diag(10))
#'
#' ## prepare input arguments for diefferent scenarios
#' lbdvec <- c(0.01,0.1,1,10,100)         # a vector of regularization parameters
#' list1 <- list(type="fixed",param=1.0)  # single regularization parameter case
#' list2 <- list(type="confidence",param=0.95) # single confidence level case
#' list3 <- list(type="BIC",param=lbdvec) # multiple regularizers with BIC selection
#'
#' ## compute with different scenarios
#' out1 <- PreEst.glasso(data, method=list1)
#' out2 <- PreEst.glasso(data, method=list2)
#' out3 <- PreEst.glasso(data, method=list3)
#'
#' ## visualize
#' par(mfrow=c(2,2), pty="s")
#' image(pracma::flipud(diag(10)),main="Original Precision")
#' image(pracma::flipud(out1$C), main="glasso::lambda=1.0")
#' image(pracma::flipud(out2$C), main="glasso::Confidence=0.95")
#' image(pracma::flipud(out3$C), main="glasso::BIC selection")
#'
#' @references
#' \insertRef{banerjee_convex_2006}{CovTools}
#'
#' \insertRef{yuan_model_2007}{CovTools}
#'
#' \insertRef{friedman_sparse_2008}{CovTools}
#'
#' @rdname PreEst.glasso
#' @export
PreEst.glasso <- function(X, method=list(type="fixed",param=1.0), parallel=FALSE){
  #-----------------------------------------------------
  ## PREPROCESSING
  #   1. typecheck : data
  fname    = "PreEst.glasso"
  checker1 = invisible_datamatrix(X, fname)
  #   2. typecheck : method
  if ((!is.list(method))||(length(method)!=2)){
    stop("* PreEst.glasso : 'method' is invalid. Please see the documentation.")
  }
  if (!(all((sort(names(method))==sort(c("type","param")))==TRUE))){
    stop("* PreEst.glasso : 'method' should have two fields; 'type' and 'param'.")
  }
  if (!(method$type %in% c("fixed","confidence","BIC"))){
    stop("* PreEst.glasso : you can opt from one of 'fixed','confidence', or 'BIC' method.")
  }
  #   3. setting for parallel computation
  if (!parallel){nCore = 1  }
  else {nCore = max(round(detectCores()/2),1)}

  #-----------------------------------------------------
  ## MAIN COMPUTATION
  #   Case 1. Fixed
  if (method$type=="fixed"){
    # parameter check
    Yuan07.lambda = method$param
    if (length(Yuan07.lambda)>1){ # vector case : it should be controlled using 'BIC'
      stop("* PreEst.glasso : 'fixed' argument only handles single regularization parameter case.")
    } else {
      if ((length(Yuan07.lambda)!=1)||(is.na(Yuan07.lambda))||(is.infinite(Yuan07.lambda))||(Yuan07.lambda<=0)){
        stop("* PreEst.glasso : for 'fixed', its correspondign value should be a positive real number.")
      }
      output = preest.Yuan07.once.matrix(X, Yuan07.lambda)
    }
  } else if (method$type=="BIC"){
    # parameter check
    Yuan07.lambda = method$param
    if (length(Yuan07.lambda)<=1){
      stop("* PreEst.glasso : 'BIC' argument handles a vector of regularization parameter case only. Use 'fixed' argument instead.")
    } else {
      if ((any(is.na(Yuan07.lambda)))||(any(is.infinite(Yuan07.lambda)))||(any(Yuan07.lambda<=0))){
        stop("* PreEst.glasso : for 'fixed', its corresponding input vector should be positive real numbers.")
      }
      output = preest.Yuan07.plgrid(X, Yuan07.lambda, nCore)
    }
  } else if (method$type=="confidence"){
    # parameter check
    Banerjee06.confidence = method$param
    if (length(Banerjee06.confidence)>1){
      stop("* PreEst.glasso : for 'confidence', it deals with a single confidence level value for each computation.")
    } else {
      if ((Banerjee06.confidence<=0)||(Banerjee06.confidence>=1)||(is.na(Banerjee06.confidence))){
        stop("* PreEst.glasso : for 'confidence', the parameter value should be in (0,1).")
      }
      output = preest.Banerjee06(X, Banerjee06.confidence)
    }
  }

  #-----------------------------------------------------
  ## RETURN OUTPUT
  return(output)
}

