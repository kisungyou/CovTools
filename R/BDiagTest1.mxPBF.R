#' One-Sample Diagonality Test by Maximum Pairwise Bayes Factor
#'
#' One-sample diagonality test can be stated with the null hypothesis
#' \deqn{H_0 : \sigma_{ij} = 0~\mathrm{for any}~i \neq j}
#' and alternative hypothesis \eqn{H_1 : ~\mathrm{not}~H_0}
#' with \eqn{\Sigma_n = (\sigma_{ij})}. Let \eqn{X_i} be the \eqn{i}-th column of data matrix. Under the maximum pairwise bayes factor framework, we have following hypothesis,
#' \deqn{H_0: a_{ij}=0\quad \mathrm{versus. } \quad  H_1: \mathrm{ not }~ H_0.}
#' The model is
#' \deqn{X_i | X_j \sim N_n( a_{ij}X_j, \tau_{ij}^2 I_n ).}
#' Under \eqn{H_0}, the prior is set as
#' \deqn{\tau_{ij}^2 \sim IG(a0, b0)} and under \eqn{H_1}, priors are
#' \deqn{ a_{ij}|\tau_{ij}^2 \sim N(0, \tau_{ij}^2/(\gamma*||X_j||^2))}
#' \deqn{\tau_{ij}^2 \sim IG(a0, b0).}
#'
#' @param data an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param a0 shape parameter for inverse-gamma prior.
#' @param b0 scale parameter for inverse-gamma prior.
#' @param gamma non-negative number. See the equation above.
#'
#' @return a named list containing: \describe{
#' \item{log.BF.mat}{ \eqn{(p\times p)} matrix of pairwise log Bayes factors.}
#' }
#'
#' @examples
#' \dontrun{
#' ## generate data from multivariate normal with trivial covariance.
#' pdim = 10
#' data = matrix(rnorm(100*pdim), nrow=100)
#'
#' ## run test
#' ## run mxPBF-based test
#' out1 = BDiagTest1.mxPBF(data)
#' out2 = BDiagTest1.mxPBF(data, a0=5.0, b0=5.0) # change some params
#'
#' ## visualize two Bayes Factor matrices
#' opar <- par(mfrow=c(1,2), pty="s")
#' image(exp(out1$log.BF.mat)[,pdim:1], main="default")
#' image(exp(out2$log.BF.mat)[,pdim:1], main="a0=b0=5.0")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{lee_maximum_2018}{CovTools}
#'
#' @export
BDiagTest1.mxPBF <- function(data, a0=2.0, b0=2.0, gamma=1.0){
  ###########################################################################
  # Preprocessing : Inputs
  # 1. data
  if (!check_datamatrix(data)){
    stop("* BDiagTest1.mxPBF : an input matrix is invalid.")
  }
  n = nrow(data)
  p = ncol(data)
  if ((nrow(data)==1)||(ncol(data)==1)){stop("* BDiagTest1.mxPBF : invalid input matrix X.")}
  # 2. extra arguments
  extra.args = list(a0=a0, b0=b0, gamma=gamma)

  ###########################################################################
  ## MAIN COMPUTATION
  output = diagtest1.mxPBF(data, extra.args)

  ## RETURN
  return(output)
}

# - Kyoungjae Lee, Lizhen Lin and David Dunson. (2018) Maximum Pairwise Bayes Factors for Covariance Structure Testing. [https://arxiv.org/abs/1809.03105]
#' @keywords internal
#' @noRd
diagtest1.mxPBF <- function(X, params){
  #################################################################
  # Params
  parnames = names(params)
  if ("a0" %in% parnames){    a0 = params$a0  } else {    a0 = 2.0} ## a0
  if ("b0" %in% parnames){    b0 = params$b0  } else {    b0 = 2.0} ## a0
  if ("gamma" %in% parnames){gamma= params$gamma}else{ gamma = 1.0}

  if ((length(a0)!=1)||(a0<=0)){stop("* mxPBF : 'a0' should be nonnegative number.")}
  if ((length(b0)!=1)||(b0<=0)){stop("* mxPBF : 'b0' should be nonnegative number.")}
  if ((length(gamma)!=1)||(gamma<=0)){stop("* mxPBF : 'gamma' should be nonnegative number.")}

  #################################################################
  ## MAIN RUN BY JAY
  p = ncol(X)
  n = nrow(X)
  log.BF.mat = matrix(0, ncol=p,nrow=p) # log Bayes factors

  for(i in 1:p){
    Xi = matrix(X[,i], ncol=1)
    for(j in (1:p)[-i]){
      Xj = matrix(X[,j], ncol=1)
      log.BF.mat[i,j] = 1/2 * log(gamma/(1+gamma)) -
        (n/2 + a0) * ( log(( sum((Xi)^2) - sum(Xi*Xj)^2/sum((Xj)^2) /(1+gamma) ) + 2*b0) - log( sum((Xi)^2) + 2*b0) )
    }
  }
  diag(log.BF.mat) = -Inf # just to fill out the diagonal parts

  output = list()
  output$log.BF.mat = log.BF.mat
  return(output)
}
