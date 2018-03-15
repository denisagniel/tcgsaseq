#' Computes gene-wise p-values from DESeq2 method using observation-wise dispersion estimates
#'
#'@param y a numeric matrix of dim \code{g x n} containing the raw RNA-seq counts for g
#'genes from \code{n} samples
#'
#'@param x a numeric design matrix of dim \code{n x p} containing the \code{p} covariates
#'to be adjusted for
#'
#'@param phi a numeric design matrix of size \code{n x K} containing the \code{K} variables
#'to be tested
#'
#'@param indiv a vector of length \code{n} containing the information for
#'attributing each sample to one of the studied individuals. Coerced
#'to be a \code{factor}
#'
#'@param ind
#'
#'@examples
#'\dontrun{
#'#rm(list=ls())
#'set.seed(123)
#'
#'##generate some fake data
#'########################
#'n <- 100
#'r <- 12
#'phi <- matrix(rep(1:3), 4, ncol=1, nrow=r)
#'sigma <- 0.4
#'b0 <- 1
#'
#'#under the null:
#'b1 <- 0
#'#under the alternative:
#'b1 <- 0.7
#'y.tilde <- b0 + b1*phi + rnorm(r, sd = sigma)
#'y <- floor(exp(t(matrix(rnorm(n*r, sd = sqrt(sigma*abs(y.tilde))), ncol=n, nrow=r) +
#'       matrix(rep(y.tilde, n), ncol=n, nrow=r))))
#'x <- matrix(1:2, ncol=1, nrow=r/2)
#'indiv=rep(1:4, each=3)
#'
#'#run test
#'temp <- deseq_fn(y, x, phi, indiv)
#'}
#'@importFrom stats as.formula
#@importFrom DESeq2 DESeqDataSetFromMatrix dispersions results estimateDispersionsGeneEst nbinomLRT
#'
#' @keywords internal
deseq_fn <- function(y, x, phi, indiv) {

  if(!requireNamespace("DESeq2", quietly=TRUE)){
    stop("Package 'DESeq2' is not available.\n  -> Try installing it from Bioconductor\n")
  }else{
    requireNamespace("S4Vectors", quietly=TRUE)

    if(is.null(colnames(phi))){
      colnames(phi) <- paste0("phi", 1:ncol(phi))
    }

    if(is.null(colnames(x))){
      colnames(x) <- paste0("x", 1:ncol(phi))
    }

    coldata_temp <- cbind.data.frame("indiv"=as.factor(indiv), phi, x)
    design_formula <- stats::as.formula(paste0("~ ", paste(colnames(x), collapse =" + "), " + ", paste(colnames(phi), collapse =" + ")))
    design_formula_reduced <- stats::as.formula(paste0("~ ", paste(colnames(x), collapse =" + ")))
    y_dsq <- DESeq2::DESeqDataSetFromMatrix(countData = y,
                                            colData = coldata_temp,
                                            design = design_formula)

    y_dsq <- DESeq2::estimateSizeFactors(y_dsq)
    y_dsq <- DESeq2::estimateDispersionsGeneEst(y_dsq)
    DESeq2::dispersions(y_dsq) <- S4Vectors::mcols(y_dsq)$dispGeneEst
    res_dsq <- DESeq2::nbinomLRT(y_dsq, reduced = design_formula_reduced)
    pvals <- DESeq2::results(res_dsq)$pvalue

    return(pvals)
  }
}
