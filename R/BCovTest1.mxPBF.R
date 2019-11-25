#' One-Sample Covariance Test using Maximum Pairwise Bayes Factor
#'
#' It performs Bayesian version of 1-sample test for Covariance where the null hypothesis is
#' \deqn{H_0 : \Sigma_n = \Sigma_0}
#' where \eqn{\Sigma_n} is the covariance of data model and \eqn{\Sigma_0} is a
#' hypothesized covariance. Denote \eqn{X_i} be the \eqn{i}-th column of data matrix.
#' Under the maximum pairwise Bayes Factor framework, we have following hypothesis,
#' \deqn{H_0: a_{ij}=0~\mathrm{ and }~\tau_{ij}=1 \quad \mathrm{versus. } \quad  H_1: \mathrm{ not }~ H_0.}
#' The model is
#' \deqn{X_i | X_j \sim N_n( a_{ij}X_j, \tau_{ij}^2 I_n )}
#' and the prior is set, under \eqn{H_1},  as
#' \deqn{ a_{ij}|\tau_{ij}^2 \sim N(0, \tau_{ij}^2/(\gamma*||X_j||^2))}
#' \deqn{\tau_{ij}^2 \sim IG(a0, b0).}
#'
#' @param data an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param Sigma0 a \eqn{(p\times p)} given covariance matrix.
#' @param a0 shape parameter for inverse-gamma prior.
#' @param b0 scale parameter for inverse-gamma prior.
#' @param gamma non-negative number. See the equation above.
#'
#' @return a named list containing: \describe{
#' \item{log.BF.mat}{a \eqn{(p\times p)} matrix of pairwise log Bayes factors.}
#' }
#'
#' @examples
#' \dontrun{
#' ## generate data from multivariate normal with trivial covariance.
#' pdim = 10
#' data = matrix(rnorm(100*pdim), nrow=100)
#'
#' ## run mxPBF-based test
#' out1 = BCovTest1.mxPBF(data)
#' out2 = BCovTest1.mxPBF(data, a0=5.0, b0=5.0) # change some params
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
BCovTest1.mxPBF <- function(data, Sigma0=diag(ncol(data)),
                            a0=2.0, b0=2.0, gamma=1.0){
  ###########################################################################
  # Preprocessing : Inputs
  # 1. data
  if (!check_datamatrix(data)){
    stop("* BCovTest1.mxPBF : an input matrix is invalid.")
  }
  n = nrow(data)
  p = ncol(data)
  if ((nrow(data)==1)||(ncol(data)==1)){stop("* BCovTest1.mxPBF : invalid input matrix X.")}
  # 2. valid Sigma0
  if ((!check_sqmat(Sigma0))||(!check_pd(Sigma0))||(!isSymmetric(Sigma0,tol=sqrt(.Machine$double.eps)))||(nrow(Sigma0)!=p)){
    stop("* BCovTest1.mxPBF : a given matrix for null hypothess 'Sigma0' is invalid.")
  }
  # 3. extra arguments
  extra.args = list(a0=a0, b0=b0, gamma=gamma)

  ###########################################################################
  # Preprocessing : Adjust the data for testing I
  scaler = get_invroot(Sigma0)
  X.centered = scale(data, center=TRUE, scale=FALSE)
  X.adjusted = (matrix(X.centered,nrow=n) %*% scaler)

  ###########################################################################
  # Main Computation
  output = bayestest1.Lee18(X.adjusted, extra.args)

  ###########################################################################
  return(output)
}


# auxiliary ---------------------------------------------------------------
# [2018.Lee] Maximum Pairwise Bayes Factors for Covariance Structure Testing
#' @keywords internal
#' @noRd
bayestest1.Lee18 <- function(X, params){
  #################################################################
  ## Params
  parnames = names(params)
  if ("a0" %in% parnames){    a0 = params$a0  } else {    a0 = 2.0} ## a0
  if ("b0" %in% parnames){    b0 = params$b0  } else {    b0 = 2.0} ## a0
  if ("gamma" %in% parnames){gamma= params$gamma}else{ gamma = 1.0}
  ## Parameter Value Warning
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
      log.BF.mat[i,j] = a0*log(b0) - lgamma(a0) +
        1/2 * log(gamma/(1+gamma)) + lgamma(n/2 + a0) +
        1/2 * sum((Xi)^2) - (n/2 + a0) * log(1/2 * ( sum((Xi)^2) - sum(Xi*Xj)^2/sum((Xj)^2) /(1+gamma) ) + b0)
    }
  }
  diag(log.BF.mat) = -Inf # just to fill out the diagonal parts
  output = list()
  output$log.BF.mat = log.BF.mat
  return(output)
}

# pack <- "CovTools"
# path <- find.package(pack)
# system(paste(shQuote(file.path(R.home("bin"), "R")),
#              "CMD", "Rd2pdf", shQuote(path)))
