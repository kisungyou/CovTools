#' A Collection of Geometric and Statistical Tools for Covariance (and Precision) Analysis
#'
#' @description
#' Covariance is of universal prevalence across various disciplines within statistics.
#' \pkg{CovTools} package aims at providing a rich collection of geometric and statistical tools
#' for a variety of inferences on covariance structures as well as its inverse called precision matrix.
#' See the sections below for a comprehensive list of functions provided from the package.
#'
#' @section Geometric Methods:
#' From inference on manifolds perspective, we have following functions,
#' \tabular{cc}{
#' \emph{name of a function} \tab \emph{description} \cr
#' \code{\link{CovDist}} \tab compute pairwise distance of covariance matrices \cr
#' \code{\link{CovMean}} \tab compute mean covariance matrix
#' }
#'
#'
#' @docType package
#' @name package-CovTools
#' @import shapes
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import Rdpack
#' @import SHT
#' @importFrom utils packageVersion
#' @importFrom pracma procrustes flipud
#' @importFrom Matrix nearPD
#' @importFrom stats cov qnorm rnorm cor qt pnorm
#' @importFrom mvtnorm rmvnorm
#' @importFrom expm expm sqrtm logm
#' @importFrom geigen geigen
#' @importFrom Rcpp evalCpp
#' @importFrom foreach "%dopar%" foreach registerDoSEQ
#' @importFrom parallel detectCores stopCluster makeCluster
#' @importFrom doParallel registerDoParallel
#' @useDynLib CovTools
NULL



