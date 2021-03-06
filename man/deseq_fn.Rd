% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deseq_fn.R
\name{deseq_fn}
\alias{deseq_fn}
\title{Computes gene-wise p-values from DESeq2 method using observation-wise dispersion estimates}
\usage{
deseq_fn(y, x, phi, indiv)
}
\arguments{
\item{y}{a numeric matrix of dim \code{g x n} containing the raw RNA-seq counts for g
genes from \code{n} samples}

\item{x}{a numeric design matrix of dim \code{n x p} containing the \code{p} covariates
to be adjusted for}

\item{phi}{a numeric design matrix of size \code{n x K} containing the \code{K} variables
to be tested}

\item{indiv}{a vector of length \code{n} containing the information for
attributing each sample to one of the studied individuals. Coerced
to be a \code{factor}}

\item{ind}{}
}
\description{
Computes gene-wise p-values from DESeq2 method using observation-wise dispersion estimates
}
\examples{
\dontrun{
#rm(list=ls())
set.seed(123)

##generate some fake data
########################
n <- 100
r <- 12
phi <- matrix(rep(1:3), 4, ncol=1, nrow=r)
sigma <- 0.4
b0 <- 1

#under the null:
b1 <- 0
#under the alternative:
b1 <- 0.7
y.tilde <- b0 + b1*phi + rnorm(r, sd = sigma)
y <- floor(exp(t(matrix(rnorm(n*r, sd = sqrt(sigma*abs(y.tilde))), ncol=n, nrow=r) +
      matrix(rep(y.tilde, n), ncol=n, nrow=r))))
x <- matrix(1:2, ncol=1, nrow=r/2)
indiv=rep(1:4, each=3)

#run test
temp <- deseq_fn(y, x, phi, indiv)
}
}
\keyword{internal}
