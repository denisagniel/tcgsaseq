#'Computes variance component test statistic for longitudinal
#'
#'This function computes an approximation of the variance component test based on the asymptotic
#'distribution of a mixture of \eqn{\chi^{2}}s using Davies method
#'from \code{\link[CompQuadForm]{davies}}
#'
#'
#'@param y a numeric matrix of dim \code{g x n} containing the raw or normalized RNAseq counts for g
#'genes from \code{n} samples.
#'
#'@param x a numeric design matrix of dim \code{n x p} containing the \code{p} covariates
#' to be adjusted for
#'
#'@param indiv a vector of length \code{n} containing the information for
#'attributing each sample to one of the studied individuals. Coerced
#'to be a \code{factor}.
#'
#'@param phi a numeric design matrix of size \code{n x K} containing the \code{K} longitudinal variables
#'to be tested (typically a vector of time points or functions of time)
#'
#'@param w a vector of length \code{n} containing the weights for the \code{n}
#'samples, corresponding to the inverse of the diagonal of the estimated covariance matrix of y.
#'
#'@param Sigma_xi a matrix of size \code{K x K} containing the covariance matrix
#'of the \code{K} random effects corresponding to \code{phi}.
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
#'
#'@seealso \code{\link[CompQuadForm]{davies}}
#'@importFrom stats cov
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
#'asymTestRes <- vc_test_asym(y, x, phi=t, w=matrix(1, ncol=ncol(y), nrow=nrow(y)),
#'                            Sigma_xi=matrix(1), indiv=rep(1:4, each=3))
#'asymTestRes$set_pval
#'
#'@importFrom CompQuadForm davies
#'@importFrom stats pchisq var
#'
#'@export
vc_test_asym <- function(y, x, indiv=rep(1,nrow(x)), phi, w, Sigma_xi = diag(ncol(phi)),
                         genewise_pvals=FALSE, homogen_traj=FALSE){


  if(homogen_traj){
    score_list <- vc_score_h(y = y, x = x, indiv = factor(indiv), phi = phi, w = w,
                             Sigma_xi = Sigma_xi)
  }else{
    score_list <- vc_score(y = y, x = x, indiv = factor(indiv), phi = phi, w = w,
                           Sigma_xi = Sigma_xi)
  }

  nindiv <- nrow(score_list$q_ext)
  ng <- ncol(score_list$q_ext)

  if(nindiv == 1){
    warning("Only 1 individual: asymptotics likely not reached - Should probably run permutation test")
    Sig_q <- matrix(1, ng, ng)
  }else{
    Sig_q <- cov(score_list$q_ext)
  }

  if (genewise_pvals) {
    if (ng == 1 & nindiv > 1) {
      gene_scores_obs <- score_list$gene_scores_unscaled/apply(score_list$q_ext, 2, stats::var)
      pv <- stats::pchisq(gene_scores_obs, df = 1, lower.tail = FALSE)
    }else if(ng > 1 & nindiv > 1){
      gene_scores_obs <- score_list$gene_scores_unscaled
      gene_lambda <- diag(Sig_q)
      pv <- unlist(mapply(FUN=CompQuadForm::davies, q=gene_scores_obs, lambda=gene_lambda, lim=15000, acc=0.0005)["Qq",])
    }else if(ng == 1 & nindiv == 1){
      gene_scores_obs <- score_list$gene_scores_unscaled
      pv <- stats::pchisq(gene_scores_obs, df = 1, lower.tail = FALSE)
    }else if(ng > 1 & nindiv == 1){
      gene_scores_obs <- score_list$gene_scores_unscaled
      pv <- stats::pchisq(gene_scores_obs, df = 1, lower.tail = FALSE)
    }else{
      stop("no gene measured/no sample included ...")
    }

    names(pv) <- rownames(y)
    ans <- list("gene_scores_obs" = gene_scores_obs, "gene_pvals" = pv)

  } else {

    lam <- try(svd(Sig_q)$d)
    if (inherits(lam, "try-error")){
      lam <- try(svd(round(Sig_q, 6))$d)
      if (inherits(lam, "try-error")){
        stop("Error in svd decomposition")
      }
    }
    dv <- CompQuadForm::davies(score_list$score, lam)
    if(dv$ifault == 1){#error
      stop("fault in the computation from CompQuadForm::davies", dv$trace)
    }
    ans <- list("set_score_obs" = score_list$score, "set_pval" = dv$Qq)
  }

  return(ans)
}

