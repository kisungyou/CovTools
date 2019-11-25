#' Two-Sample Covariance Test by Cai and Ma (2013)
#'
#' Given two sets of data, it performs 2-sample test for equality of covariance matrices where
#' the null hypothesis is
#' \deqn{H_0 : \Sigma_1 = \Sigma_2}
#' where \eqn{\Sigma_1} and \eqn{\Sigma_2} represent true (unknown) covariance
#' for each dataset based on a procedure proposed by Cai and Ma (2013).
#' If \code{statistic} \eqn{>} \code{threshold}, it rejects null hypothesis.
#'
#' @param X an \eqn{(m\times p)}  matrix where each row is an observation from the first dataset.
#' @param Y an \eqn{(n\times p)} matrix where each row is an observation from the second dataset.
#' @param alpha level of significance.
#'
#' @return a named list containing \describe{
#' \item{statistic}{a test statistic value.}
#' \item{threshold}{rejection criterion to be compared against test statistic.}
#' \item{reject}{a logical; \code{TRUE} to reject null hypothesis, \code{FALSE} otherwise.}
#' }
#'
#' @examples
#' ## generate 2 datasets from multivariate normal with identical covariance.
#' pdim  = 5
#' data1 = matrix(rnorm(100*pdim), ncol=pdim)
#' data2 = matrix(rnorm(150*pdim), ncol=pdim)
#'
#' ## run test
#' CovTest2.2013Cai(data1, data2)
#'
#' @references
#' \insertRef{cai_optimal_2013}{CovTools}
#'
#' @export
CovTest2.2013Cai <- function(X, Y, alpha=0.05){
  #---------------------------------------------------------------------------------
  ## PREPROCESSING ON INPUTS AND PARAMETERS
  # 1. valid data matrix
  if (!check_datamatrix(X)){stop("* CovTest2.2013Cai : an input data matrix X is invalid.")}
  if (!check_datamatrix(Y)){stop("* CovTest2.2013Cai : an input data matrix Y is invalid.")}
  if ((nrow(X)==1)||(ncol(X)==1)){stop("* CovTest2.2013Cai : invalid input matrix X.")}
  if ((nrow(Y)==1)||(ncol(Y)==1)){stop("* CovTest2.2013Cai : invalid input matrix Y.")}
  if (ncol(X)!=ncol(Y)){stop("* CovTest2.2013Cai : two input matrices X and Y should have same number of columns.")}
  # 2. alpha
  if ((length(alpha)!=1)||(alpha<=0)||(alpha>=1)){stop("* CovTest2 : 'alpha' should be a real number in (0,1).")}

  #---------------------------------------------------------------------------------
  ## MAIN COMPUTATION
  output = test2.Cai13(X, Y, alpha)

  #---------------------------------------------------------------------------------
  ## RETURN OUTPUT
  return(output)
}


# [2013.Cai] Two Sample Test ----------------------------------------------
# Want to test the null hypothesis H0: Sigma1 = Sigma2.
###########################################
# Input
# X : n1 X p data matrix from some distribution having covariance Sigma1
# Y : n2 X p data matrix from some distribution having covariance Sigma2
# method : Cai, Liu and Xia (2013, JASA)
# alpha : significance level
#
# Output
# This function returns a list of three arguments.
# (1) statistic : test statistic
# (2) threshold : determination criterion
# (3) reject    : result of test. TRUE means we reject H0: Sigma1 = Sigma2.
###########################################
#' @keywords internal
#' @noRd
test2.Cai13 <- function(X, Y, alpha){
  # parameter setting
  n1 = nrow(X)
  n2 = nrow(Y)
  p  = ncol(X)
  # elementary computation : sample covariance with new df
  Sigma1.hat = cov(X)*(n1-1)/n1
  Sigma2.hat = cov(Y)*(n2-1)/n2
  # mean estimation
  bar.X = colMeans(X)
  theta1.hat = matrix(0, nrow=p, ncol=p)
  for(i in 1:p){
    for(j in 1:p){
      theta1.hat[i,j] = sum( ((X[,i]-bar.X[i])*(X[,j]-bar.X[j]) - Sigma1.hat[i,j])^2 )/n1
    }
  }
  bar.Y = colMeans(Y)
  theta2.hat = matrix(0, nrow=p, ncol=p)
  for(i in 1:p){
    for(j in 1:p){
      theta2.hat[i,j] = sum( ((Y[,i]-bar.Y[i])*(Y[,j]-bar.Y[j]) - Sigma2.hat[i,j])^2 )/n2
    }
  }
  # statistic and results
  M.mat = (Sigma1.hat - Sigma2.hat)^2 / (theta1.hat/n1 + theta2.hat/n2)
  Mn = max(M.mat) # test statistic
  q.alpha  = - log(8*pi) - 2*log( log((1-alpha)^{-1}) )
  threshold= q.alpha + 4*log(p) - log(log(p))
  res = (Mn > threshold) # Tony's test
  # return results
  return(list(statistic=Mn, threshold=threshold, reject=res))
}
