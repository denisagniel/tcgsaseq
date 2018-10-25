#'Estimating the False Discovery Rate
#'
#'This function uses the permutation plug-in method to estimate the FDR.
#'
#'
#'@param gene_scores_perm a numeric matrix of size \code{G x n} containing the permuted gene-wise scores
#'for \code{G} genes with \code{n_perm} permutations.
#'
#'@param gene_scores_obs a vector of length \code{n} containing the observed gene-wise scores.
#'
#'@param n_perm the number of perturbations. Default is \code{1000}.
#'
#'@return A vector of estimating discrete false discovery rates
#'
#'@references J. Li and R. Tibshirani (2013). Finding consistent patterns: A nonparametric approach
#'for identifying differential expression in RNA-seq data, \emph{Statistical Methods in Medical Research},
#'22(5): 519-536
#'
#'@references J.D. Storey & R. Tibshirani (2003) Statistical significance for genomewide studies.
#'\emph{Proceedings of the National Academy of Sciences} 100(&6): 9440–9445.
#'
#'
#'@importFrom stats median
#'@export
dsFDR <- function(gene_scores_perm, gene_scores_obs, n_perm){

  V <- rep(0,length(gene_scores_obs))
  R <- rep(0,length(gene_scores_obs))
  FDR <- rep(0,length(gene_scores_obs))
  q <- stats::median(gene_scores_perm)
  pi_0 <- 2*mean(gene_scores_obs > q)
  for (c in 1:length(gene_scores_obs)){
    V[c] <- sum(gene_scores_perm >= gene_scores_obs[c])/n_perm
    R[c] <- sum(gene_scores_obs >= gene_scores_obs[c])
    FDR[c] <- min(1, pi_0*V[c]/R[c])
  }
  return(FDR)
}


