#'Precision weights accounting for heteroscedasticity in RNA-seq count data
#'
#'Implementation of the procedure described in Law \emph{et al.} for estimating
#'precision weights from RNA-seq data.
#'
#'@param y a matrix of size \code{G x n} containing the raw RNA-seq counts or
#'preprocessed expressions from \code{n} samples for \code{G} genes.
#'
#'@param x a matrix of size \code{n x p} containing the model covariates from
#'\code{n} samples (design matrix).
#'
#'@param preprocessed a logical flag indicating whether the expression data have
#'already been preprocessed (e.g. log2 transformed). Default is \code{FALSE}, in
#'which case \code{y} is assumed to contain raw counts and is normalized into
#'log(counts) per million.
#'
#'@param lowess_span smoother span for the lowess function, between 0 and 1.
#'This gives the proportion of points in the plot which influence the smooth at
#'each value. Larger values give more smoothness. Default is \code{0.5}.
#'
#'@param R library.size (optional, important to provide if
#'\code{preprocessed = TRUE}). Default is \code{NULL}
#'
#'@return a vector of length \code{n} containing the computed precision weights
#'
#'@seealso \code{\link{lowess}} \code{\link{approxfun}}
#'\code{\link[limma]{voom}}
#'
#'@references Law, C. W., Chen, Y., Shi, W., & Smyth, G. K. (2014). voom:
#'Precision weights unlock linear model analysis tools for RNA-seq read counts.
#'\emph{Genome Biology}, 15(2), R29.
#'
#'@examples
#'set.seed(123)
#'
#'G <- 10000
#'n <- 12
#'p <- 2
#'
#'y <- sapply(1:n, FUN=function(x){rnbinom(n=G, size=0.07, mu=200)})
#'x <- sapply(1:p, FUN=function(x){rnorm(n=n, mean=n, sd=1)})
#'
#'my_w <-  voom_weights(y, x)
#'plot_weights(my_w)
#'if (requireNamespace('limma', quietly = TRUE)) {
#'  w_voom <- limma::voom(counts=y, design=x, plot=TRUE)
#'  #slightly faster, same results
#'  all.equal(my_w$weights, w_voom$weights)
#'}
#'
#'if(interactive()){
#'#microbenchmark::microbenchmark(limma::voom(counts=t(y), design=x,
#'#                               plot=FALSE), voom_weights(x, y),
#'#                               times=30)
#'}
#'
#'@author Boris Hejblum
#'
#'@import ggplot2
#'@importFrom stats approxfun lowess
#'@export


voom_weights <- function(y, x, preprocessed = FALSE,
                         lowess_span = 0.5, R = NULL) {
    
    ## dimensions check------
    
    stopifnot(is.matrix(y))
    stopifnot(is.matrix(x))
    
    n <- ncol(y)  # the number of samples measured
    g <- nrow(y)  # the number of genes measured
    qq <- ncol(x)  # the number of covariates
    stopifnot(nrow(x) == n)
    
    # transforming rna to log-counts per million (lcpm)
    if (preprocessed) {
        y_lcpm <- t(y)
    } else {
        y_lcpm <- t(apply(y, MARGIN = 2, function(v) {
            log2((v + 0.5)/(sum(v) + 1) * 10^6)
        }))
    }
    
    # library size
    if(is.null(R)){
        R <- colSums(y, na.rm = TRUE)
    }
    
    # fitting OLS to the lcpm
    B_ols <- solve(t(x) %*% x) %*% t(x) %*% y_lcpm
    mu <- x %*% B_ols
    resi <- y_lcpm - mu
    s <- sqrt(apply(resi, MARGIN = 2, crossprod)/(nrow(x) - ncol(x)))
    s_rs <- sqrt(s)
    
    # average lcpm
    y_bar <- colMeans(y_lcpm, na.rm = TRUE)
    
    # average log-counts
    r_tilde <- y_bar + mean(log2(R + 1)) - log2(10^6)
    
    # fitted log-counts
    lambda <- mu + log2(R + 1) - log2(10^6)
    
    
    
    # lowess fit for predicted square root sd
    observed <- which(rowSums(y) != 0)  #removing genes never observed
    lowess_fit <- stats::lowess(x = r_tilde[observed], y = s_rs[observed],
                                f = lowess_span)
    f_interp <- stats::approxfun(lowess_fit, rule = 2, ties=mean)
    
    # weights
    wg <- t(apply(lambda, MARGIN = 2, f_interp)^(-4))
    
    colnames(wg) <- rownames(y_lcpm)
    
    return(list("weights" = wg,
                "plot_utilities" = list("reverse_trans" = identity,
                                        "method" = "voom",
                                        "smth" = lowess_fit,
                                        "gene_based" = NA,
                                        "mu" = r_tilde,
                                        "v" = s_rs
                )
    ))
}


