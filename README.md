
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

| function name     | authors                   | description                                              |
|-------------------|---------------------------|:---------------------------------------------------------|
| `CovEst.adaptive` | Cai and Liu (2011)        | adaptive thresholding                                    |
| `CovEst.hard`     | Bickel and Levina (2008)  | hard thresholding                                        |
| `CovEst.hardPD`   | Fan et al. (2013)         | hard thresholding under positive-definiteness constraint |
| `CovEst.nearPD`   | Qi and Sun (2006)         | nearest positive-definite matrix projection              |
| `CovEst.soft`     | Antoniadis and Fan (2001) | soft thresholding                                        |
| `CovEst.2003LW`   | Ledoit and Wolf (2003)    | linear shrinkage estimation                              |
| `CovEst.2010OAS`  | Chen et al. (2010)        | oracle approximation shrinkage                           |
| `CovEst.2010RBLW` | Chen et al. (2010)        | Rao-Blackwell Ledoit-Wolf estimation                     |

### (2) Estimation : Precision

| function name         | authors                    | description                                      |
|-----------------------|----------------------------|:-------------------------------------------------|
| `PreEst.2014An`       | An et al. (2014)           | banded precision estimation via bandwidth test   |
| `PreEst.2014Banerjee` | Banerjee and Ghosal (2014) | Bayesian estimation of a banded precision matrix |
| `PreEst.2017Lee`      | Lee and Lee (2017)         | Bayesian estimation of a banded precision matrix |
| `PreEst.glasso`       | Friedman et al. (2008)     | graphical lasso                                  |

### (3) Hypothesis Test : 1-sample

| function name             | authors                  | description                                       |
|---------------------------|--------------------------|:--------------------------------------------------|
| `BCovTeset1.mxPBF`        | Lee et al. (2018)        | Bayesian test using Maximum Pairwise Bayes Factor |
| `CovTest1.2013Cai`        | Cai and Ma (2013)        | Test by Cai and Ma                                |
| `CovTest1.2014Srivastava` | Srivastava et al. (2014) | Test by Srivastava, Yanagihara, and Kubokawa      |

### (4) Hypothesis Test : 2-sample

| function name      | authors           | description        |
|--------------------|-------------------|:-------------------|
| `CovTest2.2013Cai` | Cai and Ma (2013) | Test by Cai and Ma |

### (5) Hypothesis Test : 1-sample Diagonal

| function name       | authors              | description                                       |
|---------------------|----------------------|:--------------------------------------------------|
| `BDiagTest1.mxPBF`  | Lee et al. (2018)    | Bayesian Test using Maximum Pairwise Bayes Factor |
| `DiagTest1.2011Cai` | Cai and Jiang (2011) | Test by Cai and Jiang                             |
| `DiagTest1.2015Lan` | Lan et al. (2015)    | Test by Lan, Luo, Tsai, Wang, and Yang            |
