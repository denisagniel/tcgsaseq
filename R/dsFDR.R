#'Estimating the False Discovery Rate
#'
#'This function uses the permutation plug-in method to estimate the FDR.
#'
#'
#'@param gene_scores_perm a numeric matrix of size \code{G x n_perm} containing the permuted gene-wise scores
#'for \code{G} genes with \code{n_perm} permutations.
#'
#'@param gene_scores_obs a vector of length \code{n} containing the observed gene-wise scores.
#'
#'@param n_perm the number of perturbations. Default is \code{1000}.
#'
#'@param use_median a logical flag indicating whether the median should be used to estimate
#' the true proportion of null features. If not, we use a range of quantiles of the permuted gene-wise scores
#' and the true proportion of null features is extrapolated from the limit of a smoothed estimate using the natural cubic spline.
#' Default is \code{TRUE}. See \emph{Storey et al.} for details
#'
#'@param doPlot a logical flag indicating whether the plot of the natural cubic spline fit should be drawn.
#' Default is \code{FALSE}. Ignored if \code{use_median} is \code{TRUE}.
#'
#'@return A vector of estimating discrete false discovery rates
#'
#'@references J. Li and R. Tibshirani (2013). Finding consistent patterns: A nonparametric approach
#'for identifying differential expression in RNA-seq data, \emph{Statistical Methods in Medical Research},
#'22(5): 519-536
#'
#'@references Storey, J. D., & Tibshirani, R. (2003). Statistical significance for genome-wide studies.
#' \emph{Proceedings of the National Academy of Sciences}, 100(16), 9440-9445.
#'
#'@importFrom stats median quantile lm predict.lm
#'@importFrom graphics abline lines plot
#'@export
#'
#'@examples
#'#rm(list=ls())
#'G <- 1000
#'nperm <- 100
#'G1 <- 300
#'G0 <- G-G1
#'
#'gene_scores_perm <- matrix(rchisq(G*nperm, df=1), ncol=nperm, nrow=G)
#'gene_scores_obs <- c(rchisq(G1, df=10), rchisq(G0, df=1))
#'
#'qvals <- dsFDR(gene_scores_perm, gene_scores_obs, nperm, use_median = FALSE, doPlot = TRUE)
#'summary(qvals)
#'eFDR_5pct <- mean(qvals[-(1:G1)]<0.05)
#'eTDR_5pct <- mean(qvals[1:G1]<0.05)
#'cat("FDR:", eFDR_5pct, " TDR:", eTDR_5pct, "\n")

dsFDR <- function(gene_scores_perm, gene_scores_obs, n_perm, use_median = TRUE, doPlot = FALSE){

  # estimating pi_0 the proportion of true null hypothesis
  if(use_median){
    q <- stats::median(gene_scores_perm)
    pi_0_hat <- min(1, 2*mean(gene_scores_obs <= q))
    #cat("pi_0_hat:", pi_0_hat, "\n")
  }else{
    lambda <- seq(0.01, 0.99, by=0.01)
    q <- stats::quantile(as.vector(gene_scores_perm), probs=lambda)
    pi_0 <- pmin(1, sapply(q, function(qq){mean(gene_scores_obs <= qq)})/lambda)
    nsdf <- as.data.frame(splines::ns(c(0,lambda, 1), df = 3))
    colnames(nsdf) <- paste0("ns", 1:ncol(nsdf))
    fm <- stats::lm(y ~ ns1 + ns2 + ns3, data=cbind.data.frame("y" = pi_0, nsdf[-c(1,nrow(nsdf)), ]))
    pi_0_hat <- min(1, stats::predict.lm(fm, newdata = nsdf[c(1,nrow(nsdf)), ])[1])

    if(doPlot){
      plot(y=pi_0, x=lambda, ylim=c(0,1), xlab="lambda", col="blue")
      lines(y=fm$fitted.values, col="red", x=lambda)
      lines(pi_0_hat, x=0, type="p", pch="*", col="red")
      abline(h=pi_0_hat, col="red", lty=2)
    }
  }

  #FFDR <- pmin(1, sapply(gene_scores_obs, function(x){sum(gene_scores_perm >= x)/sum(gene_scores_obs>= x )})*pi_0_hat/nperm)
  #this is not faster - see microbenchmark::microbenchmark

  V <- rep(NA,length(gene_scores_obs))
  R <- rep(NA,length(gene_scores_obs))
  FDR <- rep(NA,length(gene_scores_obs))

  for (c in 1:length(gene_scores_obs)){
    V[c] <- sum(gene_scores_perm[c,] >= gene_scores_obs[c])
    R[c] <- sum(gene_scores_obs >= gene_scores_obs[c])
    FDR[c] <- min(1, pi_0_hat*V[c]/R[c])
  }
  return(FDR)
}
