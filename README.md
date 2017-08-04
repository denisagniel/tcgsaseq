
<!-- README.md is generated from README.Rmd. Please edit that file -->
`tcgsaseq`
==========

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/tcgsaseq)](https://cran.r-project.org/package=tcgsaseq) [![Travis-CI Build Status](https://travis-ci.org/denisagniel/tcgsaseq.svg?branch=master)](https://travis-ci.org/denisagniel/tcgsaseq) [![Downloads](https://cranlogs.r-pkg.org/badges/tcgsaseq?color=blue)](https://www.r-pkg.org/pkg/tcgsaseq)

Overview
--------

`tcgsaseq` is a package for analyzing RNA-seq data. The 2 main functions of the package are `varseq` and `tcgsa_seq`:

-   **Gene-wise Differential Analysis of RNA-seq data** can be performed using the function `varseq`.
-   **Gene Set Analysis of RNA-seq data** can be performed using the function `tcgsa_seq`.

The method implemented in this package is detailed in the following article:

> D Agniel, BP Hejblum, Variance component score test for time-course gene set analysis of longitudinal RNA-seq data, 2017, [*Biostatistics*](https://doi.org/10.1093/biostatistics/kxx005). [arXiv:1605.02351](https://arxiv.org/abs/1605.02351v4) [DOI:10.1093/biostatistics/kxx005](https://doi.org/10.1093/biostatistics/kxx005)

Installation
------------

The easiest way to get `tcgsaseq` is to install it from [CRAN](https://cran.r-project.org/package=tcgsaseq):

``` r
install.packages("tcgsaseq")
```

Or to get the development version from [GitHub](https://github.com/denisagniel/tcgsaseq):

``` r
#install.packages("devtools")
devtools::install_github("denisagniel/tcgsaseq")
```

-- Denis Agniel and Boris Hejblum
