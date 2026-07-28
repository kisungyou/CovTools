library(CovTools)

pair_values <- getFromNamespace("eigenvalues_spd_pair", "CovTools")

A <- diag(c(2, 6))
B <- diag(c(1, 2))
values <- pair_values(A, B)
stopifnot(is.double(values), !is.complex(values))
stopifnot(isTRUE(all.equal(sort(values), c(2, 3), tolerance=1e-12)))

scalar_values <- pair_values(matrix(6), matrix(2))
stopifnot(isTRUE(all.equal(scalar_values, 3, tolerance=1e-12)))

repeated_values <- pair_values(diag(rep(4, 2)), diag(rep(2, 2)))
stopifnot(isTRUE(all.equal(repeated_values, c(2, 2), tolerance=1e-12)))

invalid <- try(pair_values(diag(c(1, -1)), diag(2)), silent=TRUE)
stopifnot(inherits(invalid, "try-error"))
stopifnot(inherits(try(pair_values(diag(2), diag(3)), silent=TRUE), "try-error"))
stopifnot(inherits(try(pair_values(matrix(c(1, 2, 0, 1), 2), diag(2)),
                           silent=TRUE), "try-error"))
stopifnot(inherits(try(pair_values(diag(2), diag(c(1, 0))), silent=TRUE),
                   "try-error"))
stopifnot(inherits(try(pair_values(matrix(c(1, NA, NA, 1), 2), diag(2)),
                           silent=TRUE), "try-error"))
stopifnot(inherits(try(pair_values(diag(2)+0i, diag(2)), silent=TRUE),
                   "try-error"))

matrices <- array(0, dim=c(2, 2, 2))
matrices[, , 1] <- A
matrices[, , 2] <- B

airm <- CovDist(matrices, method="AIRM")
stopifnot(is.double(airm), !is.complex(airm))
airm_expected <- sqrt(sum(log(c(2, 3))^2))
stopifnot(isTRUE(all.equal(airm[1, 2], airm_expected, tolerance=1e-12)))

matrices_reversed <- matrices[, , 2:1, drop=FALSE]
airm_reversed <- CovDist(matrices_reversed, method="AIRM")
stopifnot(isTRUE(all.equal(airm[1, 2], airm_reversed[1, 2],
                           tolerance=1e-12)))

identical_matrices <- array(A, dim=c(2, 2, 2))
stopifnot(CovDist(identical_matrices, method="AIRM")[1, 2] < 1e-12)

kldm <- CovDist(matrices, method="KLDM")
kldm_expected <- sum((sqrt(c(2, 3))-1/sqrt(c(2, 3)))^2)/2
stopifnot(isTRUE(all.equal(kldm[1, 2], kldm_expected, tolerance=1e-12)))

set.seed(20260728)
M <- matrix(rnorm(16), nrow=4)
N <- matrix(rnorm(16), nrow=4)
A <- crossprod(M)+diag(4)
B <- crossprod(N)+diag(4)
values <- pair_values(A, B)
reference <- eigen(solve(B, A), only.values=TRUE)$values
stopifnot(isTRUE(all.equal(sort(values), sort(Re(reference)), tolerance=1e-10)))

ill_conditioned <- diag(c(1, 1e-10))
values <- pair_values(diag(c(2, 3)), ill_conditioned)
stopifnot(all(is.finite(values)), all(values>0))
