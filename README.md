
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CovTools

<!-- badges: start -->

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/CovTools?color=green)](https://cran.r-project.org/package=CovTools)
[![Travis build
status](https://travis-ci.org/kisungyou/CovTools.svg?branch=master)](https://travis-ci.org/kisungyou/CovTools)
[![](https://cranlogs.r-pkg.org/badges/CovTools)](https://cran.r-project.org/package=CovTools)
<!-- badges: end -->

Covariance is of universal prevalence across various disciplines within
statistics. This package aims at providing a rich collection of
geometric and statistical tools for a variety of inferences on
**covariance** structures as well as its inverse called **precision**
matrix. See the package help file by `help("package-CovTools")` in R
console for the list of available functions.

## Installation

You can install the released version of CovTools from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("CovTools")
```

or the development version from github:

``` r
## install.packages("devtools")
## library(devtools)
devtools::install_github("kisungyou/CovTools")
```

## List of Available Methods

We offer various methods for covariance and symmetric positive-definite
matrices. Below is the list of functions implemented in our package.

### (0) Elementary Operations

| function name | description                                                         |
|---------------|:--------------------------------------------------------------------|
| `CovDist`     | computes pairwise distance for symmetric positive-definite matrices |
| `CovMean`     | estimate mean/average covariance matrix                             |

### (1) Estimation : Covariance

| function name     | authors                                                                             | description                                              |
|-------------------|-------------------------------------------------------------------------------------|:---------------------------------------------------------|
| `CovEst.adaptive` | [Cai and Liu (2011)](https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10560) | adaptive thresholding                                    |
| `CovEst.hard`     | [Bickel and Levina (2008)](https://projecteuclid.org/euclid.aos/1231165180)         | hard thresholding                                        |
| `CovEst.hardPD`   | [Fan et al. (2013)](https://doi.org/10.1111/rssb.12016)                             | hard thresholding under positive-definiteness constraint |
| `CovEst.nearPD`   | [Qi and Sun (2006)](https://doi.org/10.1137/050624509)                              | nearest positive-definite matrix projection              |
| `CovEst.soft`     | [Antoniadis and Fan (2001)](https://doi.org/10.1198/016214501753208942)             | soft thresholding                                        |
| `CovEst.2003LW`   | [Ledoit and Wolf (2003)](https://doi.org/10.1016/S0927-5398(03)00007-0)             | linear shrinkage estimation                              |
| `CovEst.2010OAS`  | [Chen et al. (2010)](https://doi.org/10.1109/TSP.2010.2053029)                      | oracle approximation shrinkage                           |
| `CovEst.2010RBLW` | [Chen et al. (2010)](https://doi.org/10.1109/TSP.2010.2053029)                      | Rao-Blackwell Ledoit-Wolf estimation                     |

### (2) Estimation : Precision

| function name         | authors                                                                | description                                      |
|-----------------------|------------------------------------------------------------------------|:-------------------------------------------------|
| `PreEst.2014An`       | [An et al. (2014)](https://doi.org/10.1093/biomet/asu006)              | banded precision estimation via bandwidth test   |
| `PreEst.2014Banerjee` | [Banerjee and Ghosal (2014)](https://doi.org/10.1214/14-EJS945)        | Bayesian estimation of a banded precision matrix |
| `PreEst.2017Lee`      | [Lee and Lee (2017)](https://arxiv.org/abs/1707.01143)                 | Bayesian estimation of a banded precision matrix |
| `PreEst.glasso`       | [Friedman et al. (2008)](https://doi.org/10.1093/biostatistics/kxm045) | graphical lasso                                  |

### (3) Hypothesis Test : 1-sample

| function name             | authors                                                                | description                                       |
|---------------------------|------------------------------------------------------------------------|:--------------------------------------------------|
| `BCovTeset1.mxPBF`        | [Lee et al. (2018)](http://arxiv.org/abs/1809.03105)                   | Bayesian test using Maximum Pairwise Bayes Factor |
| `CovTest1.2013Cai`        | [Cai and Ma (2013)](https://doi.org/10.3150/12-BEJ455)                 | Test by Cai and Ma                                |
| `CovTest1.2014Srivastava` | [Srivastava et al. (2014)](https://doi.org/10.1016/j.jmva.2014.06.003) | Test by Srivastava, Yanagihara, and Kubokawa      |

### (4) Hypothesis Test : 2-sample

| function name      | authors                                                | description        |
|--------------------|--------------------------------------------------------|:-------------------|
| `CovTest2.2013Cai` | [Cai and Ma (2013)](https://doi.org/10.3150/12-BEJ455) | Test by Cai and Ma |

### (5) Hypothesis Test : 1-sample Diagonal

| function name       | authors                                                                              | description                                       |
|---------------------|--------------------------------------------------------------------------------------|:--------------------------------------------------|
| `BDiagTest1.mxPBF`  | [Lee et al. (2018)](http://arxiv.org/abs/1809.03105)                                 | Bayesian Test using Maximum Pairwise Bayes Factor |
| `DiagTest1.2011Cai` | [Cai and Jiang (2011)](http://projecteuclid.org/euclid.aos/1305292044)               | Test by Cai and Jiang                             |
| `DiagTest1.2015Lan` | [Lan et al. (2015)](http://www.tandfonline.com/doi/abs/10.1080/07350015.2014.923317) | Test by Lan, Luo, Tsai, Wang, and Yang            |
