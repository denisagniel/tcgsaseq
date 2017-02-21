#'Permutation-based variance component test statistic for longitudinal RNA-seq data
#'
#'This function computes an approximation of the Variance Component test for a
#'mixture of \eqn{\chi^{2}}s using permutations. This is preferable to the
#'asymptotic approximation for small sample sizes
#'
#'
#'@param y a numeric matrix of dim \code{g x n} containing the raw RNAseq counts for g
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
#'@param n_perm the number of perturbations
#'
#'@param genewise_pvals a logical flag indicating whether genewise pvalues should be returned. Default
#'is \code{FALSE} in which case geneset p-value is computed and returned instead.
#'
#'@param homogen_traj a logical flag indicating whether trajectories should be considered homogeneous.
#'Default is \code{FALSE} in which case trajectories are not only tested for trend, but also for heterogeneity.
#'
#'@return A list with the following elements when the set p-value is computed :\itemize{
#'   \item \code{set_score_obs}: the approximation of the observed set score
#'   \item \code{set_pval}: the associated set p-value
#' }
#' or a list with the following elements when genewise pvalues are computed:\itemize{
#'   \item \code{gene_scores_obs}: vector of approximating the observed genewise scores
#'   \item \code{gene_pvals}: vector of associated genewise p-values
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
vc_test_perm <- function(y, x, indiv=rep(1,nrow(x)), phi, w, Sigma_xi = diag(ncol(phi)),
                         n_perm=1000, genewise_pvals=FALSE, homogen_traj=FALSE){

  n_samples <- ncol(y)
  if(is.null(colnames(y))){
    colnames(y) <- 1:n_samples
  }
  n_genes <- nrow(y)

  indiv_fact <- factor(indiv)

  if(homogen_traj){
    vc_score_2use <- vc_score_h
  }else{
    vc_score_2use <- vc_score
  }

  strat_sampling <- function(fact){
    res <- numeric(length(fact))
    for(l in levels(fact)){
      original_index <- which(fact==l)
      res[original_index] <- as.numeric(sample(as.character(original_index)))
    }
    return(res)
  }

  score_list_obs <- vc_score_2use(y = y, x = x, indiv = indiv_fact, phi = phi, w = w, Sigma_xi = Sigma_xi)

  if(genewise_pvals){
    gene_scores_obs <- score_list_obs$gene_scores_unscaled

    gene_scores_perm <- list()
    for(b in 1:n_perm){
      ## permute samples within indiv
      ## tried to pertub genes by gene (to avoid inducing dependency between genes, doesn't change a thing)
      # gene_scores_perm[[b]] <- sapply(1:n_genes, function(i){
      #   perm_index <- strat_sampling(indiv_fact)
      #   res <- vc_score_2use(y = y[i, perm_index, drop=FALSE], x = x[perm_index, , drop=FALSE],
      #                        indiv = indiv_fact, phi = phi, w = w[i, perm_index, drop = FALSE],
      #                        Sigma_xi = Sigma_xi)$gene_scores_unscaled
      #   return(res)})

      perm_index <- strat_sampling(indiv_fact)
      gene_scores_perm[[b]] <- vc_score_2use(y = y[, perm_index, drop=FALSE], x = x[perm_index, , drop=FALSE],
                             indiv = indiv_fact, phi = phi, w = w[, perm_index, drop = FALSE],
                             Sigma_xi = Sigma_xi)$gene_scores_unscaled
    }
    #browser()
    pvals <- 1-rowMeans(sapply(gene_scores_perm, function(x){x<gene_scores_obs}))
    pvals2 <- 1-rowMeans(do.call(cbind, gene_scores_perm)<gene_scores_obs)
    #hist(pvals)
    ans <- list("gene_scores_obs" = gene_scores_obs, "gene_pvals" = pvals)
  }else{
    scores_perm <- numeric(n_perm)
    for(b in 1:n_perm){
      ## permute samples within indiv
      perm_index <- strat_sampling(indiv_fact)
      scores_perm[b] <- vc_score_2use(y = y[, perm_index, drop=FALSE], x = x[perm_index, , drop=FALSE],
                                      indiv = indiv_fact, phi = phi, w = w[, perm_index, drop = FALSE], Sigma_xi = Sigma_xi)$score
    }
    pval <- 1-sum(scores_perm < score_list_obs$score)/n_perm
    ans <- list("set_score_obs" = score_list_obs$score, "set_pval" = pval)
  }

}
