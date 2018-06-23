**tranSurv**
------------

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![minimal R version](https://img.shields.io/badge/R%3E%3D-3.4.4-6666ff.svg)](https://cran.r-project.org/) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/tranSurv)](https://cran.r-project.org/package=tranSurv) [![packageversion](https://img.shields.io/badge/Package%20version-1.1.6-orange.svg?style=flat-square)](commits/master) [![Last-changedate](https://img.shields.io/badge/last%20change-2018--06--23-yellowgreen.svg)](/commits/master)

<!-- README.md is generated from README.Rmd. Please edit that file -->
#### Transformation model for left-truncated right-censored survival data

***tranSurv*** implements methods for survival analysis under a dependent truncation and independent right censoring via a structural transformation method. The package is still under active development. \*\*\*

A package that estimates survival curve under a dependent truncation and independent right censoring via a structural transformation method. The package also includes hypothesis test of quasi-independence based on the conditional Kendall's tau of Martin and Betensky (2005) and two versions of the inverse probability weighted Kendall's tau of Austin and Betensky (2014).

Installation
------------

Install and load the package from CRAN using

``` r
install.packages("tranSurv")
library(tranSurv)
```

Install and load the package from GitHub using

``` r
devtools::install_github("stc04003/tranSurv")
library(tranSurv)
```
