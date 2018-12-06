#'Estimating the False Discovery Rate
#'
#'This function uses the permutation plug-in method to estimate the FDR.
#'It requires scaled-test statistics so that permutation are comparable from one test to another.
#'
#'
#'@param gene_scores_perm a numeric matrix of size \code{G x n_perm} containing the permuted gene-wise scores
#'for \code{G} genes with \code{n_perm} permutations.
#'
#'@param gene_scores_obs a vector of length \code{n} containing the observed gene-wise scores.
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
#'@importFrom graphics abline lines plot legend
#'@export
#'
#'@examples
#'\dontrun{
#'#rm(list=ls())
#'G <- 1000
#'nperm <- 500
#'G1 <- 0.3*G
#'G0 <- G-G1
#'gene_scores_perm <- matrix(rchisq(G*nperm, df=1), ncol=nperm, nrow=G)
#'gene_scores_obs <- c(rchisq(G1, df=1), rchisq(G0, df=1))
#'
#'qvals <- dsFDR(gene_scores_perm, gene_scores_obs, use_median = FALSE, doPlot = TRUE)
#'summary(qvals)
#'
#'qvals <- qvals[!is.na(qvals)]
#'eFDR_5pct <- sum(qvals[-(1:G1)]<0.05)/sum(qvals < 0.05)
#'eTDR_5pct <- sum(qvals[1:G1]<0.05)/sum(qvals < 0.05)
#'cat("FDR:", eFDR_5pct, " TDR:", eTDR_5pct, "\n")
#'plot(y = sapply(seq(0, 1, by=0.001), function(x){sum(qvals[-(1:G1)] < x)/sum(qvals < x)}),
#'     x = seq(0, 1, by=0.001),
#'     type = "l", xlab = "Nominal FDR level", ylab = "Empirical FDR", col = "red", lwd = 2,
#'     ylim = c(0,1))
#'abline(a = 0, b = 1, lty = 2)
#'
#'res <- list()
#'G <- 1000
#'nperm <- 1000
#'G1 <- 0.3*G
#'G0 <- G-G1
#'for(i in 1:100){
#'cat(i, "/100\n", sep="")
#'gene_scores_obs <- c(rchisq(G1, df=1), rchisq(G0, df=1))
#'gene_scores_perm <- matrix(rchisq(G*nperm, df=1), ncol=nperm, nrow=G)
#'qvals <- dsFDR(gene_scores_perm, gene_scores_obs, use_median = TRUE, doPlot = FALSE)
#'res[[i]] <- sapply(seq(0, 1, by=0.01), function(x){sum(qvals < x)})
#'}
#'}

dsFDR <- function(gene_scores_perm, gene_scores_obs, use_median = TRUE, doPlot = FALSE){

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
      legend("bottomright", c("estimate", "spline fit"), lty=c(2,1), col="red", pch=c("*", ""))
    }
  }


  #V_li <- sapply(gene_scores_obs, function(threshold){sum(rowMeans(gene_scores_perm > threshold))})
  #R_li <- sapply(gene_scores_obs, function(threshold){sum(gene_scores_obs > threshold)})
  #FDR_li <- pmin(1, pi_0_hat*V_li/R_li)
  FDR_li <- numeric(length(gene_scores_obs))
  for(i in 1:length(gene_scores_obs)){
    FDR_li[i] <- min(1, pi_0_hat*sum(rowMeans(gene_scores_perm > gene_scores_obs[i]))/sum(gene_scores_obs > gene_scores_obs[i]))
    message(paste0(i, "/", length(gene_scores_obs), "\n"))
  }

  #V_scaled <- rowMeans(gene_scores_perm >= gene_scores_obs)
  #R_scaled <- sapply(gene_scores_obs, function(x){mean(gene_scores_obs >= x)})
  #FDR <- pmin(1, pi_0_hat*V_scaled/R_scaled)


  return(FDR_li)
}
