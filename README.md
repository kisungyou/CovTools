
<!-- README.md is generated from README.Rmd. Please edit that file -->
CovTools
========

<!-- badges: start -->
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/CovTools?color=green)](https://cran.r-project.org/package=CovTools) [![Travis build status](https://travis-ci.org/kyoustat/CovTools.svg?branch=master)](https://travis-ci.org/kyoustat/CovTools) [![](https://cranlogs.r-pkg.org/badges/CovTools)](https://cran.r-project.org/package=CovTools) <!-- badges: end -->

Covariance is of universal prevalence across various disciplines within statistics. This package aims at providing a rich collection of geometric and statistical tools for a variety of inferences on **covariance** structures as well as its inverse called **precision** matrix. See the package help file by `help("package-CovTools")` in R console for the list of available functions.

Installation
------------

You can install the released version of CovTools from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("CovTools")
```

or the development version from github:

``` r
## install.packages("devtools")
## library(devtools)
devtools::install_github("kyoustat/CovTools")
```

List of Available Methods
-------------------------

We offer various methods for covariance and symmetric positive-definite matrices. Below is the list of functions implemented in our package.

### (0). Elementary Operations

<table style="width:51%;">
<colgroup>
<col width="22%" />
<col width="29%" />
</colgroup>
<thead>
<tr class="header">
<th>function name</th>
<th align="left">description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>CovDist</code></td>
<td align="left">computes pairwise distance for symmetric positive-definite matrices</td>
</tr>
<tr class="even">
<td><code>CovMean</code></td>
<td align="left">estimate mean/average covariance matrix</td>
</tr>
</tbody>
</table>

### (1). Estimation : Covariance

<table style="width:78%;">
<colgroup>
<col width="22%" />
<col width="26%" />
<col width="29%" />
</colgroup>
<thead>
<tr class="header">
<th>function name</th>
<th>authors</th>
<th align="left">description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>CovEst.adaptive</code></td>
<td><a href="https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10560">Cai and Liu (2011)</a></td>
<td align="left">adaptive thresholding</td>
</tr>
<tr class="even">
<td><code>CovEst.hard</code></td>
<td><a href="https://projecteuclid.org/euclid.aos/1231165180">Bickel and Levina (2008)</a></td>
<td align="left">hard thresholding</td>
</tr>
<tr class="odd">
<td><code>CovEst.hardPD</code></td>
<td><a href="https://doi.org/10.1111/rssb.12016">Fan et al. (2013)</a></td>
<td align="left">hard thresholding under positive-definiteness constraint</td>
</tr>
<tr class="even">
<td><code>CovEst.nearPD</code></td>
<td><a href="https://doi.org/10.1137/050624509">Qi and Sun (2006)</a></td>
<td align="left">nearest positive-definite matrix projection</td>
</tr>
<tr class="odd">
<td><code>CovEst.soft</code></td>
<td><a href="https://doi.org/10.1198/016214501753208942">Antoniadis and Fan (2001)</a></td>
<td align="left">soft thresholding</td>
</tr>
<tr class="even">
<td><code>CovEst.2003LW</code></td>
<td><a href="https://doi.org/10.1016/S0927-5398(03)00007-0">Ledoit and Wolf (2003)</a></td>
<td align="left">linear shrinkage estimation</td>
</tr>
<tr class="odd">
<td><code>CovEst.2010OAS</code></td>
<td><a href="https://doi.org/10.1109/TSP.2010.2053029">Chen et al. (2010)</a></td>
<td align="left">oracle approximation shrinkage</td>
</tr>
<tr class="even">
<td><code>CovEst.2010RBLW</code></td>
<td><a href="https://doi.org/10.1109/TSP.2010.2053029">Chen et al. (2010)</a></td>
<td align="left">Rao-Blackwell Ledoit-Wolf estimation</td>
</tr>
</tbody>
</table>

### (2). Estimation : Precision

<table style="width:78%;">
<colgroup>
<col width="22%" />
<col width="26%" />
<col width="29%" />
</colgroup>
<thead>
<tr class="header">
<th>function name</th>
<th>authors</th>
<th align="left">description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>PreEst.2014An</code></td>
<td><a href="https://doi.org/10.1093/biomet/asu006">An et al. (2014)</a></td>
<td align="left">banded precision estimation via bandwidth test</td>
</tr>
<tr class="even">
<td><code>PreEst.2014Banerjee</code></td>
<td><a href="https://doi.org/10.1214/14-EJS945">Banerjee and Ghosal (2014)</a></td>
<td align="left">Bayesian estimation of a banded precision matrix</td>
</tr>
<tr class="odd">
<td><code>PreEst.2017Lee</code></td>
<td><a href="https://arxiv.org/abs/1707.01143">Lee and Lee (2017)</a></td>
<td align="left">Bayesian estimation of a banded precision matrix</td>
</tr>
<tr class="even">
<td><code>PreEst.glasso</code></td>
<td><a href="https://doi.org/10.1093/biostatistics/kxm045">Friedman et al. (2008)</a></td>
<td align="left">graphical lasso</td>
</tr>
</tbody>
</table>

### (3). Hypothesis Test : 1-sample

<table style="width:78%;">
<colgroup>
<col width="22%" />
<col width="26%" />
<col width="29%" />
</colgroup>
<thead>
<tr class="header">
<th>function name</th>
<th>authors</th>
<th align="left">description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>BCovTeset1.mxPBF</code></td>
<td><a href="http://arxiv.org/abs/1809.03105">Lee et al. (2018)</a></td>
<td align="left">Bayesian test using Maximum Pairwise Bayes Factor</td>
</tr>
<tr class="even">
<td><code>CovTest1.2013Cai</code></td>
<td><a href="https://doi.org/10.3150/12-BEJ455">Cai and Ma (2013)</a></td>
<td align="left">Test by Cai and Ma</td>
</tr>
<tr class="odd">
<td><code>CovTest1.2014Srivastava</code></td>
<td><a href="https://doi.org/10.1016/j.jmva.2014.06.003">Srivastava et al. (2014)</a></td>
<td align="left">Test by Srivastava, Yanagihara, and Kubokawa</td>
</tr>
</tbody>
</table>

### (4). Hypothesis Test : 2-sample

| function name      | authors                                                | description        |
|--------------------|--------------------------------------------------------|:-------------------|
| `CovTest2.2013Cai` | [Cai and Ma (2013)](https://doi.org/10.3150/12-BEJ455) | Test by Cai and Ma |

### (5). Hypothesis Test : 1-sample Diagonal

<table style="width:78%;">
<colgroup>
<col width="22%" />
<col width="26%" />
<col width="29%" />
</colgroup>
<thead>
<tr class="header">
<th>function name</th>
<th>authors</th>
<th align="left">description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>BDiagTest1.mxPBF</code></td>
<td><a href="http://arxiv.org/abs/1809.03105">Lee et al. (2018)</a></td>
<td align="left">Bayesian Test using Maximum Pairwise Bayes Factor</td>
</tr>
<tr class="even">
<td><code>DiagTest1.2011Cai</code></td>
<td><a href="http://projecteuclid.org/euclid.aos/1305292044">Cai and Jiang (2011)</a></td>
<td align="left">Test by Cai and Jiang</td>
</tr>
<tr class="odd">
<td><code>DiagTest1.2015Lan</code></td>
<td><a href="http://www.tandfonline.com/doi/abs/10.1080/07350015.2014.923317">Lan et al. (2015)</a></td>
<td align="left">Test by Lan, Luo, Tsai, Wang, and Yang</td>
</tr>
</tbody>
</table>
