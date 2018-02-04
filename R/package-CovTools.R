#' A Collection of Geometric and Statistical Tools for Covariance Analysis
#'
#' Covariance is of universal prevalence across various disciplines within statistics.
#' The \pkg{CovTools} aims at providing a rich collection of geometric and inferential tools
#' for convenient analysis of covariance structures. Following is a list of functions,
#' \tabular{cc}{
#' \emph{name of a function} \tab \emph{description} \cr
#' \code{\link{CovDist}} \tab compute pairwise distance of covariance matrices \cr
#' \code{\link{CovEst}}  \tab estimate covariance matrix \cr
#' \code{\link{CovEst.auto}}\tab estimate covariance matrix with automatic parameter tuning \cr
#' \code{\link{CovMean}}\tab compute mean covariance matrix \cr
#' \code{\link{CovTest1}} \tab 1-sample tests for covariance matrix \cr
#' \code{\link{CovTest2}} \tab 2-sample tests for covariance matrices \cr
#' \code{\link{PreEst}} \tab estimate an inverse covariance matrix \cr
#' \code{\link{PreEst.auto}}\tab estimate an inverse covariance matrix with automatic parameter tuning
#' }
#'
#' @docType package
#' @name CovTools-package
#' @import shapes
#' @import parallel
#' @import doParallel
#' @import foreach
#' @importFrom pracma procrustes flipud
#' @importFrom Matrix nearPD
#' @importFrom stats cov qnorm rnorm cor qt
#' @importFrom mvtnorm rmvnorm
#' @importFrom expm expm sqrtm logm
#' @importFrom geigen geigen
#' @importFrom Rcpp evalCpp
#' @importFrom foreach "%dopar%" foreach registerDoSEQ
#' @importFrom parallel detectCores stopCluster makeCluster
#' @importFrom doParallel registerDoParallel
#' @useDynLib CovTools
NULL



