library(CovTools)

set.seed(20260728)
X <- matrix(rnorm(60*4), nrow=60, ncol=4)
location.shift <- c(7, -11, 2.5, 19)
X.shifted <- sweep(X, 2L, location.shift, FUN="+", check.margin=FALSE)

cai.original <- DiagTest1.2011Cai(X)
cai.shifted <- DiagTest1.2011Cai(X.shifted)
stopifnot(identical(names(cai.original), c("statistic", "threshold", "reject")))
stopifnot(is.logical(cai.original$reject), length(cai.original$reject)==1L)
stopifnot(isTRUE(all.equal(cai.original, cai.shifted, tolerance=1e-10)))

lan.original <- DiagTest1.2015Lan(X)
lan.shifted <- DiagTest1.2015Lan(X.shifted)
stopifnot(isTRUE(all.equal(lan.original, lan.shifted, tolerance=1e-10)))

constant.column <- X
constant.column[, 1L] <- 1
lan.constant.error <- try(DiagTest1.2015Lan(constant.column), silent=TRUE)
stopifnot(inherits(lan.constant.error, "try-error"))
stopifnot(grepl("positive sample variation", as.character(lan.constant.error),
                fixed=TRUE))

bayes.original <- BDiagTest1.mxPBF(X)
bayes.shifted <- BDiagTest1.mxPBF(X.shifted)
stopifnot(isTRUE(all.equal(
  bayes.original$log.BF.mat,
  bayes.shifted$log.BF.mat,
  tolerance=1e-10
)))

too.small <- matrix(rnorm(3*3), nrow=3, ncol=3)
sample.size.error <- try(CovTest1.2014Srivastava(too.small), silent=TRUE)
stopifnot(inherits(sample.size.error, "try-error"))
stopifnot(grepl("at least 4 observations", as.character(sample.size.error),
                fixed=TRUE))

minimum.size <- matrix(rnorm(4*3), nrow=4, ncol=3)
minimum.size.result <- CovTest1.2014Srivastava(minimum.size)
stopifnot(is.finite(minimum.size.result$statistic))
stopifnot(is.logical(minimum.size.result$reject))

cv.fun <- getFromNamespace("thr1.multiple", "CovTools")
cv.environment <- new.env(parent=environment(cv.fun))
cv.environment$makeCluster <- function(...){
  stop("makeCluster must not be called by sequential cross-validation.")
}
environment(cv.fun) <- cv.environment

hard.from.cov <- function(S, thr){
  output <- S
  diagonal <- diag(output)
  output[abs(S)<=thr] <- 0
  diag(output) <- diagonal
  output
}

backend.before <- foreach::getDoParName()
set.seed(314159)
cv.original <- cv.fun(
  X=X,
  nCV=4L,
  nCore=1L,
  func_S=hard.from.cov,
  thrvec=c(0.05, 0.10, 0.20)
)
backend.after <- foreach::getDoParName()
stopifnot(identical(backend.before, backend.after))
stopifnot(nrow(cv.original$CV)==3L)
stopifnot(all(is.finite(cv.original$CV$CVscore)))

set.seed(314159)
cv.repeated <- cv.fun(
  X=X,
  nCV=4L,
  nCore=1L,
  func_S=hard.from.cov,
  thrvec=c(0.05, 0.10, 0.20)
)
stopifnot(isTRUE(all.equal(cv.original, cv.repeated, tolerance=0)))

invalid.split <- try(
  cv.fun(
    X=matrix(rnorm(4*3), nrow=4, ncol=3),
    nCV=2L,
    nCore=1L,
    func_S=hard.from.cov,
    thrvec=c(0.05, 0.10)
  ),
  silent=TRUE
)
stopifnot(inherits(invalid.split, "try-error"))
stopifnot(grepl("at least two observations", as.character(invalid.split),
                fixed=TRUE))
