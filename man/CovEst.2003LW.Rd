% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CovEst.2003LW.R
\name{CovEst.2003LW}
\alias{CovEst.2003LW}
\title{Covariance Estimation with Linear Shrinkage}
\usage{
CovEst.2003LW(X, target = NULL)
}
\arguments{
\item{X}{an \eqn{(n\times p)} matrix where each row is an observation.}

\item{target}{target matrix \eqn{F}. If \code{target=NULL}, \emph{constant correlation model} estimator is used. If \code{target} is specified as a qualified matrix, it is used instead.}
}
\value{
a named list containing: \describe{
\item{S}{a \eqn{(p\times p)} covariance matrix estimate.}
\item{delta}{an estimate for convex combination weight according to the relevant theory.}
}
}
\description{
Ledoit and Wolf (2003, 2004) proposed a linear shrinkage strategy to estimate covariance matrix
with an application to portfolio optimization. An optimal covariance is written as a convex combination as follows,
\deqn{\hat{\Sigma} = \delta \hat{F} + (1-\delta) \hat{S}}
where \eqn{\delta \in (0,1)} a control parameter/weight, \eqn{\hat{S}} an empirical covariance matrix, and \eqn{\hat{F}} a \emph{target} matrix.
Although authors used \eqn{F} a highly structured estimator, we also enabled an arbitrary target matrix to be used as long as it's symmetric
and positive definite of corresponding size.
}
\examples{
## CRAN-purpose small computation
# set a seed for reproducibility
set.seed(11)

#  small data with identity covariance
pdim      <- 5
dat.small <- matrix(rnorm(20*pdim), ncol=pdim)

#  run the code with highly structured estimator
out.small <- CovEst.2003LW(dat.small)

#  visualize
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,3), pty="s")
image(diag(5)[,pdim:1], main="true cov")
image(cov(dat.small)[,pdim:1], main="sample cov")
image(out.small$S[,pdim:1], main="estimated cov")
par(opar)

\dontrun{
## want to see how delta is determined according to
#  the number of observations we have.
nsamples = seq(from=5, to=200, by=5)
nnsample = length(nsamples)

#  we will record two values; delta and norm difference
vec.delta = rep(0, nnsample)
vec.normd = rep(0, nnsample)
for (i in 1:nnsample){
  dat.norun <- matrix(rnorm(nsamples[i]*pdim), ncol=pdim) # sample in R^5
  out.norun <- CovEst.2003LW(dat.norun)                   # run with default

  vec.delta[i] = out.norun$delta
  vec.normd[i] = norm(out.norun$S - diag(pdim),"f")       # Frobenius norm
}

# let's visualize the results
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(nsamples, vec.delta, lwd=2, type="b", col="red", main="estimated deltas")
plot(nsamples, vec.normd, lwd=2, type="b", col="blue",main="Frobenius error")
par(opar)
}

}
\references{
\insertRef{ledoit_improved_2003}{CovTools}

\insertRef{ledoit_well-conditioned_2004}{CovTools}

\insertRef{ledoit_honey_2004}{CovTools}
}
