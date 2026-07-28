library(CovTools)

is_spd <- function(A){
  is.matrix(A) &&
    isTRUE(isSymmetric(A, tol=1e-10)) &&
    all(is.finite(A)) &&
    min(eigen(A, symmetric=TRUE, only.values=TRUE)$values) > 0
}

expect_error <- function(expr){
  error <- tryCatch({
    force(expr)
    NULL
  }, error=identity)
  stopifnot(inherits(error, "error"))
}

# Banerjee (2014): a selected bandwidth smaller than upperK must still fill
# the complete precision matrix, and tied modes must resolve to one scalar.
set.seed(101)
X_banerjee <- matrix(rnorm(60*6), nrow=60, ncol=6)
first_mode_prior <- function(k){
  if (k == 1) 0 else -1e6
}
banerjee_stein <- PreEst.2014Banerjee(
  X_banerjee, upperK=3, logpi=first_mode_prior, loss="Stein"
)
banerjee_squared <- PreEst.2014Banerjee(
  X_banerjee, upperK=3, logpi=first_mode_prior, loss="Squared"
)
stopifnot(
  is_spd(banerjee_stein$C),
  is_spd(banerjee_squared$C),
  all(diag(banerjee_stein$C) > 0),
  all(diag(banerjee_squared$C) > 0)
)

banerjee_logjg <- getFromNamespace("preest.Banerjee14.logJG", "CovTools")
X_banerjee_centered <- scale(X_banerjee, center=TRUE, scale=FALSE)
banerjee_offsets <- vapply(
  1:3,
  function(k) as.numeric(banerjee_logjg(k, X_banerjee_centered, 10)),
  numeric(1)
)
tied_banerjee_prior <- function(k){
  -banerjee_offsets[k]
}
banerjee_tied <- PreEst.2014Banerjee(
  X_banerjee, upperK=3, logpi=tied_banerjee_prior, loss="Stein"
)
banerjee_k1 <- PreEst.2014Banerjee(
  X_banerjee, upperK=1, logpi=function(k) 0, loss="Stein"
)
stopifnot(isTRUE(all.equal(banerjee_tied$C, banerjee_k1$C, tolerance=1e-10)))
expect_error(PreEst.2014Banerjee(X_banerjee, upperK=1.5))

# An et al. (2014): alpha affects selection, Holm handles no rejection, and
# centering makes the precision estimator invariant to a location shift.
set.seed(2)
X_an <- matrix(rnorm(80*6), nrow=80, ncol=6)
an_low_alpha <- PreEst.2014An(
  X_an, upperK=3, algorithm="Bonferroni", alpha=0.01
)
an_high_alpha <- PreEst.2014An(
  X_an, upperK=3, algorithm="Bonferroni", alpha=0.9
)
stopifnot(an_low_alpha$optk == 1L, an_high_alpha$optk == 3L)
an_shifted <- PreEst.2014An(
  X_an+10, upperK=3, algorithm="Bonferroni", alpha=0.9
)
stopifnot(
  is_spd(an_high_alpha$C),
  isTRUE(all.equal(an_shifted$C, an_high_alpha$C, tolerance=1e-10)),
  identical(an_shifted$optk, an_high_alpha$optk)
)

set.seed(5)
X_an_holm <- matrix(rnorm(100*5), nrow=100, ncol=5)
an_holm <- PreEst.2014An(
  X_an_holm, upperK=2, algorithm="Holm", alpha=0.01
)
stopifnot(an_holm$optk == 1L, is_spd(an_holm$C))
expect_error(PreEst.2014An(X_an, upperK=1.5))
expect_error(PreEst.2014An(matrix(rnorm(8*10), 8, 10), upperK=4))

# Lee (2017): tied modes are scalar, output remains SPD, and invalid
# bandwidths are rejected before a singular regression is attempted.
set.seed(303)
X_lee <- matrix(rnorm(50*6), nrow=50, ncol=6)
lee_logpost <- getFromNamespace("Lee17.logpostpi", "CovTools")
X_lee_centered <- scale(X_lee, center=TRUE, scale=FALSE)
lee_offsets <- vapply(
  1:3,
  function(k) lee_logpost(k, X_lee_centered, function(k) 0),
  numeric(1)
)
tied_lee_prior <- function(k){
  -lee_offsets[k]
}
lee_tied <- PreEst.2017Lee(X_lee, upperK=3, logpi=tied_lee_prior)
lee_k1 <- PreEst.2017Lee(X_lee, upperK=1, logpi=function(k) 0)
stopifnot(
  is_spd(lee_tied$C),
  isTRUE(all.equal(lee_tied$C, lee_k1$C, tolerance=1e-10))
)
expect_error(PreEst.2017Lee(X_lee, upperK=1.5))
expect_error(PreEst.2017Lee(matrix(rnorm(8*10), 8, 10), upperK=6))

# Graphical lasso: return an exactly sparse SPD matrix, compute the stated
# per-observation BIC, and run the sequential BIC path without a socket.
set.seed(404)
X_glasso <- matrix(rnorm(80*6), nrow=80, ncol=6)
glasso_fixed <- PreEst.glasso(
  X_glasso, method=list(type="fixed", param=1)
)
off_diagonal <- glasso_fixed$C[row(glasso_fixed$C) != col(glasso_fixed$C)]
stopifnot(is_spd(glasso_fixed$C), any(off_diagonal == 0))

lambda <- 0.25
glasso_score <- getFromNamespace("preest.Yuan07.once.BIC", "CovTools")(
  X_glasso, lambda
)
glasso_lambda <- PreEst.glasso(
  X_glasso, method=list(type="fixed", param=lambda)
)$C
sample_covariance <- cov(X_glasso)
logdet <- determinant(glasso_lambda, logarithm=TRUE)
degrees_freedom <- sum(
  abs(glasso_lambda[upper.tri(glasso_lambda, diag=TRUE)]) >
    sqrt(.Machine$double.eps)
)
expected_score <- -as.numeric(logdet$modulus) +
  sum(glasso_lambda*sample_covariance) +
  log(nrow(X_glasso))*degrees_freedom/nrow(X_glasso)
stopifnot(isTRUE(all.equal(glasso_score, expected_score, tolerance=1e-10)))

glasso_bic <- PreEst.glasso(
  X_glasso,
  method=list(type="BIC", param=c(0.1, 0.5, 1)),
  parallel=FALSE
)
stopifnot(
  is_spd(glasso_bic$C),
  identical(glasso_bic$BIC$lambda, c(0.1, 0.5, 1)),
  all(is.finite(glasso_bic$BIC$BIC))
)

parallel_capable <- tryCatch({
  probe_cluster <- parallel::makeCluster(2)
  parallel::stopCluster(probe_cluster)
  TRUE
}, error=function(e) FALSE)
if (parallel_capable){
  backend_name <- foreach::getDoParName()
  backend_workers <- foreach::getDoParWorkers()
  glasso_parallel <- getFromNamespace("preest.Yuan07.plgrid", "CovTools")(
    X_glasso, c(0.1, 0.5), nCore=2
  )
  stopifnot(
    is_spd(glasso_parallel$C),
    identical(foreach::getDoParName(), backend_name),
    identical(foreach::getDoParWorkers(), backend_workers)
  )
}
