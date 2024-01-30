
<!-- README.md is generated from README.Rmd. Please edit that file -->

# interep

> **Inter**action Analysis of **Rep**eated Measure Data

<!-- badges: start -->

<!-- [![CRAN](https://www.r-pkg.org/badges/version/interep)](https://cran.r-project.org/package=interep) -->

<!-- [![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/interep)](http://www.r-pkg.org/pkg/interep) -->

[![Travis build
status](https://travis-ci.org/feizhoustat/interep.svg?branch=master)](https://app.travis-ci.com/feizhoustat/interep)
[![CRAN
status](https://www.r-pkg.org/badges/version/interep)](https://CRAN.R-project.org/package=interep)
[![Codecov test
coverage](https://codecov.io/gh/feizhoustat/interep/branch/master/graph/badge.svg)](https://app.codecov.io/gh/feizhoustat/interep?branch=master)
<!-- badges: end -->

Extensive penalized variable selection methods have been developed in the past two decades for analyzing high dimensional omics data, such as gene expressions, single nucleotide polymorphisms (SNPs), copy number variations (CNVs) and others. However, lipidomics data have been rarely investigated by using high dimensional variable selection methods. This package incorporates our recently developed penalization procedures to conduct interaction analysis for high dimensional lipidomics data with repeated measurements. The core module of this package is developed in C++. The development of this software package and the associated statistical methods have been partially supported by an Innovative Research Award from Johnson Cancer Research Center, Kansas State University.

## How to install

  - Released versions of interep are available on CRAN
    [(link)](https://cran.r-project.org/package=interep), and can be
    installed within R via

<!-- end list -->

    install.packages("interep")

## Example

    library(interep)
    data("dat")
    ## Load the environment factors, lipid factors and the response
    e=dat$e
    g=dat$z
    y=dat$y
    ## Initial value for the coefficient vector
    beta0=dat$coef
    ## True nonzero coefficients
    index=dat$index
    b = interep(e, g, y,beta0,corre="e",pmethod="mixed",lam1=dat$lam1, lam2=dat$lam2,maxits=30)
    ## Cut off the noise
    b[abs(b)<0.05]=0
    ## Compute TP and FP
    pos = which(b != 0)
    tp = length(intersect(index, pos))
    fp = length(pos) - tp
    list(tp=tp, fp=fp)



## Methods

This package provides implementation for methods proposed in

  - Zhou, F., Ren,J., Li, G., Jiang, Y., Li, X., Wang, W. and Wu, C. (2019). Penalized Variable Selection for Lipid–        Environment Interactions in a Longitudinal Lipidomics Study. Genes. 10(12), 1002
