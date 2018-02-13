#' A Collection of Geometric and Statistical Tools for Covariance Analysis
#'
#' Covariance is of universal prevalence across various disciplines within statistics.
#' The \pkg{CovTools} aims at providing a rich collection of geometric and statistical tools
#' for a variety of inferences on covariance structures. Following is a list of functions,
#' \tabular{cc}{
#' \emph{name of a function} \tab \emph{description} \cr
#' \code{\link{CovDist}} \tab compute pairwise distance of covariance matrices \cr
#' \code{\link{CovMean}}\tab compute mean covariance matrix \cr
#' \code{\link{CovTest1}} \tab 1-sample tests for covariance matrix \cr
#' \code{\link{CovTest2}} \tab 2-sample tests for covariance matrices
#' }
#'
#' Also, below is a list of functions for \emph{estimating} covariance matrices,
#' \tabular{cc}{
#' \emph{name of a function}\tab \emph{description} \cr
#' \code{\link{CovEst.adaptive}} \tab Adaptive Thresholding \cr
#' \code{\link{CovEst.hard}} \tab Hard Thresholding \cr
#' \code{\link{CovEst.hardPD}} \tab Hard Thresholding under Positive-Definiteness Constraint \cr
#' \code{\link{CovEst.nearPD}} \tab Nearest Positive-Definite Matrix Projection \cr
#' \code{\link{CovEst.soft}} \tab Soft Thresholding
#' }
#' and precision - inverse covariance - matrices,
#' \tabular{cc}{
#' \emph{name of a function}\tab \emph{description} \cr
#' \code{\link{PreEst.glasso}} \tab Graphical Lasso
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



