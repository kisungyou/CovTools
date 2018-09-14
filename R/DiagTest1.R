#' One-Sample Diagonality Tests for Covariance Matrix
#'
#' Given data, \code{DiagTest1} performs Frequentist and Bayesian version of 1-sample test for diagonal entries of a Covariance matrix where
#' the null hypothesis is
#' \deqn{H_0 : \sigma_{ij} = 0~\mathrm{for any}~i \neq j}
#' and alternative hypothesis is \eqn{H_1 : ~\mathrm{not}~H_0}
#' with \eqn{\Sigma_n = (\sigma_{ij})}.
#'
#' @param data an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param alpha level of significance.
#' @param method a name of test.
#' @param ... extra arguments to be passed along for each procedure. See below for details.
#' \tabular{lll}{
#' \emph{parameter} \tab \emph{method} \tab \emph{description} \cr
#' \code{a0} \tab \code{"mxPBF"} \tab hyperparameter (see below for details) \cr
#' \code{b0} \tab \code{"mxPBF"} \tab hyperparameter (see below for details) \cr
#' \code{gamma} \tab \code{"mxPBF"} \tab hyperparameter (see below for details)
#' }
#'
#' @return a named list containing one of followings,
#' \tabular{lll}{
#' \emph{element} \tab \emph{method} \tab \emph{description} \cr
#' \code{statistic} \tab Frequentist \tab a test statistic value. \cr
#' \code{threshold} \tab Frequentist \tab rejection criterion to be compared against test statistic. \cr
#' \code{reject} \tab Frequentist \tab a logical; \code{TRUE} to reject null hypothesis, \code{FALSE} otherwise. \cr
#' \code{log.BF.mat} \tab \code{"mxPBF"} \tab a \eqn{(p\times p)} matrix of pairwise log Bayes factors.
#' }
#'
#' @section mxPBF:
#' Let \eqn{X_i} be the \eqn{i}-th column of data matrix. Under the maximum pairwise bayes factor framework, we have following hypothesis,
#' \deqn{H_0: a_{ij}=0\quad \mathrm{versus. } \quad  H_1: \mathrm{ not }~ H_0.}
#' The model is
#' \deqn{X_i | X_j \sim N_n( a_{ij}X_j, \tau_{ij}^2 I_n ).}
#' Under \eqn{H_0}, the prior is set as
#' \deqn{\tau_{ij}^2 \sim IG(a0, b0)} and under \eqn{H_1}, priors are
#' \deqn{ a_{ij}|\tau_{ij}^2 \sim N(0, \tau_{ij}^2/(\gamma*||X_j||^2))}
#' \deqn{\tau_{ij}^2 \sim IG(a0, b0).}
#'
#' @examples
#' \dontrun{
#' ## generate data from multivariate normal with trivial covariance.
#' data = matrix(rnorm(100*5), nrow=100)
#'
#' ## run test
#' DiagTest1(data, method="Cai11")
#' DiagTest1(data, method="Lan15")
#' DiagTest1(data, method="mxPBF")
#' DiagTest1(data, method="mxPBF", a0=5.0, b0=5.0) # change hyperparameters for mxPBF
#' }
#'
#' @references
#' \insertRef{cai_limiting_2011}{CovTools}
#'
#' \insertRef{lan_testing_2015}{CovTools}
#'
#' \insertRef{lee_maximum_2018}{CovTools}
#'
#' @export
DiagTest1 <- function(data, alpha=0.05, method=c("Lan15","Cai11","mxPBF"), ...){
  ###########################################################################
  # Preprocessing : Inputs
  # 1. data
  if (!check_datamatrix(data)){
    stop("* DiagTest1 : an input matrix is invalid.")
  }
  n = nrow(data)
  p = ncol(data)
  if ((nrow(data)==1)||(ncol(data)==1)){stop("* BayesTest1 : invalid input matrix X.")}
  # 2. method
  method = match.arg(method)
  # 3. alpha
  if ((length(alpha)!=1)||(alpha<=0)||(alpha>=1)){
    stop("* DiagTest1 : 'alpha' should be a real number in (0,1).")
  }
  # 4. extra arguments
  extra.args = (list(...))

  ###########################################################################
  ## MAIN COMPUTATION
  output = switch(method,
                  Lan15 = diagtest1.Lan15(data, alpha),
                  Cai11 = diagtest1.Cai11(data, alpha),
                  mxPBF = diagtest1.mxPBF(data, extra.args)
                  )

  ###########################################################################
  ## RETURN OUTPUT
  return(output)
}
