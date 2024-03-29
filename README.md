
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `tcgsaseq`

<!-- badges: start -->

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/tcgsaseq)](https://cran.r-project.org/package=tcgsaseq)
[![R-CMD-check](https://github.com/denisagniel/tcgsaseq/workflows/R-CMD-check/badge.svg)](https://github.com/denisagniel/tcgsaseq/actions)
[![Coveralls test
coverage](https://coveralls.io/repos/github/denisagniel/tcgsaseq/badge.svg)](https://coveralls.io/github/denisagniel/tcgsaseq)
[![Downloads](https://cranlogs.r-pkg.org/badges/tcgsaseq?color=blue)](https://www.r-pkg.org/pkg/tcgsaseq)
<!-- badges: end -->

<img src="man/figures/lifecycle-deprecated.svg" width="20%" />

## WARNING: deprecation

`tcgsaseq` capabilities and features have been incorporated into the
[`Bioconductor` package
`dearseq`](https://bioconductor.org/packages/release/bioc/html/dearseq.html).

`tcgsaseq` is the legacy companion package of the 2017 Agniel & Hejblum
[*Biostatistics*](https://academic.oup.com/biostatistics/article-abstract/18/4/589/3065599)
article. It has been archived on CRAN and will not be further developped
or extended.

If you are interested in using tcgsaseq, you should probably rather look
into the newer and more comprehensive `dearseq` package.

## Overview

`tcgsaseq` is a package for analyzing RNA-seq data. The 2 main functions
of the package are `varcompseq` and `tcgsa_seq`:

-   **Gene-wise Differential Analysis of RNA-seq data** can be performed
    using the function `varcompseq`.
-   **Gene Set Analysis of RNA-seq data** can be performed using the
    function `tcgsa_seq`.

The method implemented in this package is detailed in the following
article:

> Agniel D & Hejblum BP (2017). Variance component score test for
> time-course gene set analysis of longitudinal RNA-seq data, 2017,
> [*Biostatistics*](https://academic.oup.com/biostatistics/article-abstract/18/4/589/3065599),
> 18(4):589-604. [arXiv:1605.02351](https://arxiv.org/abs/1605.02351v4)
> [DOI:
> 10.1093/biostatistics/kxx005](https://doi.org/10.1093/biostatistics/kxx005)

## Installation

The easiest way to get `tcgsaseq` is to install it from
[GitHub](https://github.com/denisagniel/tcgsaseq):

``` r
#install.packages("devtools")
devtools::install_github("denisagniel/tcgsaseq")
```

– Denis Agniel and Boris Hejblum
