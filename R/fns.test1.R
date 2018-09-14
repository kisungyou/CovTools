# [2013.Cai] Optimal hypothesis testing for high-dimensional covar --------
# Tony Cai and Zongming Ma. (2013) Optimal hypothesis testing for high dimensional covariance matrices.
# If (Tn > Trej), reject Null Hypothesis (Sigma=I)
#' @keywords internal
#' @noRd
test1.Cai13 <- function(X, alpha){
  # 3. parameters and hXiXj
  n = nrow(X)
  p = ncol(X)
  hXiXj = rcpptest1_cai11(X)
  # 4. post-adjusting
  statistic = (hXiXj*2)/(n*(n-1)) # test statistic
  threshold = (qnorm(1-alpha))*2*sqrt((p*(p+1))/(n*(n-1)))
  reject = (statistic > threshold)
  # 5. return
  return(list(statistic=statistic, threshold=threshold, reject=reject))
}

# [2014.Srivastava] Tests for covariance matrices in high dimension with less sample size.
#' @keywords internal
#' @noRd
test1.Srivastava14 <- function(X, alpha){
  Sigma0 = diag(ncol(X))
  p = ncol(X)
  n = nrow(X)
  z.val = qnorm(1 - alpha)
  bar.X = matrix(colMeans(X), nrow = p, ncol = 1)
  Sn = cov(X)
  Yn = X - matrix(c(bar.X), nrow=n, ncol=p, byrow=T)
  a1.hat = sum(diag(Sn)) / p
  a2.hat = ( (n-1)*(n-2)*sum(diag( (n-1)^2*Sn%*%Sn )) - n*(n-1)*sum((rowSums(Yn^2))^2) + (sum(diag( (n-1)*Sn )))^2 ) /( p*n*(n-1)*(n-2)*(n-3) )
  T2 = (n-1)/2 * ( a2.hat - 2*a1.hat + 1 )

  res = list()
  res$statistic = T2
  res$threshold = z.val
  res$reject = ( T2 > z.val )
  return( res )
}

