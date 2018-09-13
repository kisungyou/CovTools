# [2013.Cai] Optimal hypothesis testing for high-dimensional covar --------
# Tony Cai and Zongming Ma. (2013) Optimal hypothesis testing for high dimensional covariance matrices.
# If (Tn > Trej), reject Null Hypothesis (Sigma=I)
#' @keywords internal
#' @noRd
test1.Cai13 <- function(X, Sigma0, alpha){
  # 1. setup
  scaler = get_invroot(Sigma0)
  # 2. mean centering
  X.centered = scale(X, center=TRUE, scale=FALSE)
  X.centered = matrix(X.centered,nrow(X)) %*% scaler
  # 3. parameters and hXiXj
  n = nrow(X)
  p = ncol(X)
  hXiXj = rcpptest1_cai11(X.centered)
  # 4. post-adjusting
  statistic = (hXiXj*2)/(n*(n-1)) # test statistic
  threshold = (qnorm(1-alpha))*2*sqrt((p*(p+1))/(n*(n-1)))
  reject = (statistic > threshold)
  # 5. return
  return(list(statistic=statistic, threshold=threshold, reject=reject))
}
