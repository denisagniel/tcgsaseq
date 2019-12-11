#'dearseq: a package to perform Differential Expression Analysis of RNA-seq data
#'
#'Analysis of RNA-seq data with variance component score
#'test accounting for data heteroscedasticity through precision weights.
#'Performs gene-wise analysis as well as gene set analysis, including for
#'complex experimental designs such as longitudinal data.
#'
#'\tabular{ll}{
#'Package: \tab dearseq\cr
#'Type: \tab Package\cr
#'Version: \tab 1.0.0\cr
#'Date: \tab 2019-12-11\cr
#'License:\tab \href{http://www.gnu.org/licenses/gpl-2.0.txt}{GPL-2}\cr
#'}
#'The two main functions of the \code{dearseq} package are
#'\code{\link{dear_seq}} and \code{\link{dgsa_seq}}.
#'
#'
#'@author Boris P. Hejblum, Denis Agniel
#'--- Maintainer: Boris P. Hejblum
#'
#'@references Agniel D & Hejblum BP (2017). Variance component score test for
#'time-course gene set analysis of longitudinal RNA-seq data,
#'\emph{Biostatistics}, 18(4):589-604.
#'\href{https://doi.org/10.1093/biostatistics/kxx005}{10.1093/biostatistics/kxx005}.
#'\href{https://arxiv.org/abs/1605.02351}{arXiv:1605.02351}.
#'
#'Gauthier M, Agniel D, Thiébaut R & Hejblum BP (2019). dearseq: a variance
#'component score test for RNA-Seq differential analysis that effectively
#'controls the false discovery rate, *bioRxiv* 635714.
#'[DOI: 10.1101/635714v1](https://www.biorxiv.org/content/10.1101/635714v1).
#'
#'@aliases dearseq-package
#'
#'@docType package
#'@name dearseq
#'
NULL
