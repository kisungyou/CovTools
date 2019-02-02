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
#' @section Statistical Methods:
#' We provide statistical methods for
#' (1) \strong{Covariance Matrix Estimation},
#' \tabular{ll}{
#' \emph{name of a function}\tab \emph{description} \cr
#' \code{\link{CovEst.adaptive}} \tab Adaptive Thresholding. \cr
#' \code{\link{CovEst.hard}} \tab Hard Thresholding. \cr
#' \code{\link{CovEst.hardPD}} \tab Hard Thresholding under Positive-Definiteness Constraint. \cr
#' \code{\link{CovEst.nearPD}} \tab Nearest Positive-Definite Matrix Projection. \cr
#' \code{\link{CovEst.soft}} \tab Soft Thresholding.
#' }
#' (2) \strong{Precision Matrix Estimation}
#' \tabular{ll}{
#' \emph{name of a function}\tab \emph{description} \cr
#' \code{\link{PreEst.2014An}} \tab Banded Precision Matrix Estimation via Bandwidth Test. \cr
#' \code{\link{PreEst.2014Banerjee}} \tab Bayesian Estimation of a Banded Precision Matrix (Banerjee 2014). \cr
#' \code{\link{PreEst.2017Lee}} \tab Bayesian Estimation of a Banded Precision Matrix (Lee 2017). \cr
#' \code{\link{PreEst.glasso}} \tab Graphical Lasso.
#' }
#' (3) \strong{1-Sample Covariance Tests}
#' \tabular{ll}{
#' \emph{name of a function} \tab \emph{description} \cr
#' \code{\link{BCovTest1.mxPBF}} \tab Bayesian Test using Maximum Pairwise Bayes Factor. \cr
#' \code{\link{CovTest1.2013Cai}} \tab Test by Cai and Ma (2013). \cr
#' \code{\link{CovTest1.2014Srivastava}} \tab Test by Srivastava, Yanagihara, and Kubokawa (2014).
#' }
#' (4) \strong{2-Sample Covariance Tests}
#' \tabular{ll}{
#' \emph{name of a function} \tab \emph{description} \cr
#' \code{\link{CovTest2.2013Cai}} \tab Test by Cai and Ma (2013).
#' }
#' (5) \strong{1-Sample Diagonality Tests}
#' \tabular{ll}{
#' \emph{name of a function} \tab \emph{description} \cr
#' \code{\link{BDiagTest1.mxPBF}} \tab Bayesian Test using Maximum Pairwise Bayes Factor. \cr
#' \code{\link{DiagTest1.2011Cai}} \tab Test by Cai and Jiang (2011). \cr
#' \code{\link{DiagTest1.2015Lan}} \tab Test by Lan et al. (2015).
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



