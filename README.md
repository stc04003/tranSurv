**tranSurv**
------------

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.4.0-6666ff.svg)](https://cran.r-project.org/)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/tranSurv)](https://cran.r-project.org/package=tranSurv)
[![packageversion](https://img.shields.io/badge/Package%20version-1.2.2-orange.svg?style=flat-square)](commits/master)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/stc04003/tranSurv?branch=master&svg=true)](https://ci.appveyor.com/project/stc04003/tranSurv)
[![Travis-CI Build
Status](https://travis-ci.org/stc04003/tranSurv.svg?branch=master)](https://travis-ci.org/stc04003/tranSurv)
[![Last-changedate](https://img.shields.io/badge/last%20change-2020--12--21-yellowgreen.svg)](/commits/master)

<!-- README.md is generated from README.Rmd. Please edit that file -->

#### Transformation model for left-truncated right-censored survival data

***tranSurv*** implements methods for survival analysis under a
dependent truncation and independent right censoring via a structural
transformation method. The package is still under active development.

------------------------------------------------------------------------

A package that estimates survival curve under a dependent truncation and
independent right censoring via a structural transformation method. The
package also includes hypothesis test of quasi-independence based on the
conditional Kendall’s tau of Martin and Betensky (2005) and two versions
of the inverse probability weighted Kendall’s tau of Austin and Betensky
(2014).

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

### Online documentation

[Online document](https://www.sychiou.com/tranSurv/index.html) includes:

-   Package vignette on [testing
    quasi-independence](https://www.sychiou.com/tranSurv/articles/tranSurv-quasi-independence.html)
-   Package vignette on [survival curve estimation with transformation
    models](https://www.sychiou.com/tranSurv/articles/tranSurv-test.html)
-   Package vignette on regression with with transformation models
    (coming up)
