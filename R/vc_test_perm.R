#'Permutation-based variance component test statistic
#'
#'This function computes an approximation of the Variance Component test for a
#'mixture of \eqn{\chi^{2}}s using permutations. This is preferable to the
#'asymptotic approximation for small sample sizes. We rely on exact p-values following
#'Phipson and Smyth, 2010 (see References).
#'
#'
#'@param y a numeric matrix of dim \code{G x n} containing the raw RNA-seq counts for G
#'genes from \code{n} samples.
#'
#'@param x a numeric design matrix of dim \code{n x p} containing the \code{p} covariates
#' to be adjusted for
#'
#'@param indiv a vector of length \code{n} containing the information for
#'attributing each sample to one of the studied individuals. Coerced
#'to be a \code{factor}.
#'
#'@param phi a numeric design matrix of size \code{n x K} containing the \code{K} variables
#'to be tested
#'
#'@param w a vector of length \code{n} containing the weights for the \code{n}
#'samples.
#'
#'@param Sigma_xi a matrix of size \code{K x K} containing the covariance matrix
#'of the \code{K} random effects.
#'
#'@param n_perm the number of perturbations. Default is \code{1000}.
#'
#'@param genewise_pvals a logical flag indicating whether gene-wise p-values should be returned. Default
#'is \code{FALSE} in which case gene-set p-value is computed and returned instead.
#'
#'@param homogen_traj a logical flag indicating whether trajectories should be considered homogeneous.
#'Default is \code{FALSE} in which case trajectories are not only tested for trend, but also for heterogeneity.
#'
#'@param na.rm logical: should missing values (including \code{NA} and \code{NaN})
#'be omitted from the calculations? Default is \code{FALSE}.
#'
#'@references Phipson B, and Smyth GK (2010). Permutation p-values should never be zero:
#'calculating exact p-values when permutations are randomly drawn. \emph{Statistical Applications
#'in Genetics and Molecular Biology}, Volume 9, Issue 1, Article 39.
#'\url{http://www.statsci.org/smyth/pubs/PermPValuesPreprint.pdf}
#'
#'@return A list with the following elements when the set p-value is computed :\itemize{
#'   \item \code{set_score_obs}: the approximation of the observed set score
#'   \item \code{set_pval}: the associated set p-value
#' }
#' or a list with the following elements when gene-wise p-values are computed:\itemize{
#'   \item \code{gene_scores_obs}: vector of approximating the observed gene-wise scores
#'   \item \code{gene_pvals}: vector of associated gene-wise p-values
#'   \item \code{ds_fdr}: vector of associated gene-wise discrete false discovery rates
#' }
#'
#'@seealso \code{\link[CompQuadForm]{davies}}
#'
#'@examples
#'#rm(list=ls())
#'set.seed(123)
#'
#'##generate some fake data
#'########################
#'n <- 100
#'r <- 12
#'t <- matrix(rep(1:3), 4, ncol=1, nrow=r)
#'sigma <- 0.4
#'b0 <- 1
#'
#'#under the null:
#'b1 <- 0
#'#under the alternative:
#'b1 <- 0.5
#'y.tilde <- b0 + b1*t + rnorm(r, sd = sigma)
#'y <- t(matrix(rnorm(n*r, sd = sqrt(sigma*abs(y.tilde))), ncol=n, nrow=r) +
#'       matrix(rep(y.tilde, n), ncol=n, nrow=r))
#'x <- matrix(1, ncol=1, nrow=r)
#'
#'#run test
#'permTestRes <- vc_test_perm(y, x, phi=t, w=matrix(1, ncol=ncol(y), nrow=nrow(y)),
#'                            indiv=rep(1:4, each=3), n_perm=50) #1000)
#'permTestRes$set_pval
#'
#'@importFrom CompQuadForm davies
#'
#'@export
vc_test_perm <- function(y, x, indiv = rep(1,nrow(x)), phi, w, Sigma_xi = diag(ncol(phi)),
                         n_perm = 1000, genewise_pvals = FALSE, homogen_traj = FALSE,
                         na.rm = FALSE){

  n_samples <- ncol(y)
  if(is.null(colnames(y))){
    colnames(y) <- 1:n_samples
  }
  indiv_fact <- factor(indiv)

  if(homogen_traj){
    vc_score_2use <- vc_score_h_perm
  }else{
    vc_score_2use <- vc_score_perm
  }

  if(is.null(indiv)){
    options(warn = -1)
    N_possible_perms <- factorial(ncol(y))
    options(warn = 0)
  }else{
    options(warn = -1)
    N_possible_perms <- prod(sapply(table(indiv), factorial))
    options(warn = 0)
  }


  score_list_res <- vc_score_2use(y = y, x = x, indiv = indiv_fact, phi = phi, w = w,
                                  Sigma_xi = Sigma_xi, na_rm = na.rm, n_perm = n_perm)

  if(genewise_pvals){
    gene_scores_obs <- score_list_res$gene_scores_unscaled
    gene_scores_perm <- score_list_res$gene_scores_unscaled_perm

    nprem_supobs <- rowSums(gene_scores_perm >= gene_scores_obs)

    #pvals_naive <- nprem_supobs/n_perm
    #pvals_u <- (nprem_supobs + 1)/(n_perm +1)
    pvals_e <- perm_pe(nprem_supobs, nperm_eff = n_perm, total_possible_nperm = N_possible_perms)
    ans <- list("gene_scores_obs" = gene_scores_obs, "gene_pvals" = pvals_e)
  }else{
    pvals_u <- (sum(score_list_res$scores_perm >= score_list_res$score) + 1)/(n_perm + 1)
    ans <- list("set_score_obs" = score_list_res$score, "set_pval" = pvals_u)

  }
}
