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

