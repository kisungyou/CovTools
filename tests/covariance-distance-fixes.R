library(CovTools)

pd_tolerance <- sqrt(.Machine$double.eps)
minimum_eigenvalue <- function(A){
  min(eigen(A, symmetric=TRUE, only.values=TRUE)$values)
}
expect_error <- function(expr){
  inherits(try(force(expr), silent=TRUE), "try-error")
}

# CovEst.hardPD: the reported high-dimensional case must return a matrix
# that passes the same positive-definiteness tolerance used by the search.
set.seed(1)
high_dimensional <- matrix(rnorm(10*20), nrow=10, ncol=20)
hard_pd <- CovEst.hardPD(high_dimensional)
stopifnot(identical(names(hard_pd), c("S", "optC")))
stopifnot(is.numeric(hard_pd$optC), length(hard_pd$optC)==1)
stopifnot(is.finite(hard_pd$optC), hard_pd$optC>=0)
stopifnot(isSymmetric(hard_pd$S, tol=1e-12))
stopifnot(minimum_eigenvalue(hard_pd$S)>=pd_tolerance)

# The returned estimate must be the hard-thresholded covariance at optC.
sample_covariance <- cov(high_dimensional)
expected_hard_pd <- sample_covariance
off_diagonal <- (row(expected_hard_pd)!=col(expected_hard_pd))
expected_hard_pd[off_diagonal &
                   (abs(expected_hard_pd)<=hard_pd$optC)] <- 0
expected_hard_pd <- (expected_hard_pd+t(expected_hard_pd))/2
stopifnot(isTRUE(all.equal(hard_pd$S, expected_hard_pd,
                           tolerance=1e-12)))

# The already-PD path must still return optC, with zero meaning no threshold.
pd_data <- rbind(c(-1, 0), c(1, 0), c(0, -1), c(0, 1))
unthresholded <- suppressMessages(CovEst.hardPD(pd_data))
stopifnot(identical(names(unthresholded), c("S", "optC")))
stopifnot(identical(unthresholded$optC, 0))
stopifnot(minimum_eigenvalue(unthresholded$S)>=pd_tolerance)

# A zero-variance column makes the requested guarantee impossible because
# hard thresholding preserves the covariance diagonal.
zero_variance <- cbind(1:5, rep(1, 5))
stopifnot(expect_error(CovEst.hardPD(zero_variance)))
stopifnot(expect_error(CovEst.hardPD(matrix(1:3, nrow=1))))

# PowerEuclidean distances must be returned as real, nonnegative matrices
# for either sign of a finite non-zero power.
A <- diag(c(4, 9))
B <- diag(c(1, 16))
power_input <- array(0, dim=c(2, 2, 2))
power_input[, , 1] <- A
power_input[, , 2] <- B

for (power in c(-2, -1, 0.5, 2)){
  distance <- CovDist(power_input, method="PowerEuclidean", power=power)
  stopifnot(is.matrix(distance), is.double(distance), !is.complex(distance))
  stopifnot(all(is.finite(distance)), all(distance>=0))
  stopifnot(isTRUE(all.equal(distance, t(distance), tolerance=1e-14)))
  stopifnot(all(diag(distance)==0))
}

negative_power_expected <- sqrt(sum((diag(A)^(-1)-diag(B)^(-1))^2))
negative_power_distance <- CovDist(power_input,
                                   method="PowerEuclidean", power=-1)
stopifnot(isTRUE(all.equal(negative_power_distance[1, 2],
                           negative_power_expected, tolerance=1e-12)))
stopifnot(expect_error(CovDist(power_input,
                              method="PowerEuclidean", power=0)))
stopifnot(expect_error(CovDist(power_input,
                              method="PowerEuclidean", power=Inf)))

# A one-dimensional covariance catches diag(scalar) dimension mistakes.
scalar_input <- array(c(4, 9), dim=c(1, 1, 2))
scalar_distance <- CovDist(scalar_input,
                           method="PowerEuclidean", power=-1)
stopifnot(identical(dim(scalar_distance), c(2L, 2L)))
stopifnot(isTRUE(all.equal(scalar_distance[1, 2],
                           abs((1/4)-(1/9)), tolerance=1e-12)))

# Every pairwise method has the same well-defined one-slice result.
one_slice <- array(diag(c(2, 3)), dim=c(2, 2, 1))
distance_methods <- c("AIRM", "Bhattacharyya", "Cholesky", "Euclidean",
                      "Hellinger", "JBLD", "KLDM", "LERM",
                      "Procrustes.SS", "Procrustes.Full",
                      "PowerEuclidean", "RootEuclidean")
for (method in distance_methods){
  distance <- if (method=="PowerEuclidean"){
    CovDist(one_slice, method=method, power=-1)
  } else {
    CovDist(one_slice, method=method)
  }
  stopifnot(identical(dim(distance), c(1L, 1L)))
  stopifnot(identical(as.double(distance), 0))
}

# Procrustes.SS must compare Cholesky factors, not the covariance matrices.
procrustes_input <- array(0, dim=c(2, 2, 2))
procrustes_input[, , 1] <- matrix(c(4, 1.2, 1.2, 2), nrow=2)
procrustes_input[, , 2] <- matrix(c(2, 0.3, 0.3, 1.5), nrow=2)
procrustes_distance <- CovDist(procrustes_input, method="Procrustes.SS")
factor_distance <- pracma::procrustes(
  t(chol(procrustes_input[, , 1])),
  t(chol(procrustes_input[, , 2]))
)$d
raw_distance <- pracma::procrustes(
  procrustes_input[, , 1],
  procrustes_input[, , 2]
)$d
stopifnot(isTRUE(all.equal(procrustes_distance[1, 2],
                           factor_distance, tolerance=1e-12)))
stopifnot(abs(procrustes_distance[1, 2]-raw_distance)>1e-6)

# Negative-power means are computed spectrally and remain finite and SPD.
negative_power_mean <- CovMean(power_input,
                               method="PowerEuclidean", power=-1)
expected_mean_diagonal <- 1/rowMeans(cbind(1/diag(A), 1/diag(B)))
stopifnot(all(is.finite(negative_power_mean)))
stopifnot(isSymmetric(negative_power_mean, tol=1e-12))
stopifnot(minimum_eigenvalue(negative_power_mean)>0)
stopifnot(isTRUE(all.equal(diag(negative_power_mean),
                           expected_mean_diagonal, tolerance=1e-12)))

scalar_mean <- CovMean(scalar_input, method="PowerEuclidean", power=-1)
stopifnot(identical(dim(scalar_mean), c(1L, 1L)))
stopifnot(isTRUE(all.equal(as.double(scalar_mean),
                           1/mean(c(1/4, 1/9)), tolerance=1e-12)))
stopifnot(expect_error(CovMean(power_input,
                              method="PowerEuclidean", power=Inf)))
