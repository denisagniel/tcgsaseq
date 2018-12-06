#'Computes variance component test statistic for longitudinal
#'
#'This function computes an approximation of the variance component test based on the asymptotic
#'distribution of a mixture of \eqn{\chi^{2}}s using Davies method
#'from \code{\link[CompQuadForm]{davies}}
#'
#'
#'@param y a numeric matrix of dim \code{g x n} containing the raw or normalized RNA-seq counts for g
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
#'@param genewise_pvals a logical flag indicating whether gene-wise p-values should be returned. Default
#'is \code{FALSE} in which case gene set p-value is computed and returned instead.
#'
#'@param homogen_traj a logical flag indicating whether trajectories should be considered homogeneous.
#'Default is \code{FALSE} in which case trajectories are not only tested for trend, but also for heterogeneity.
#'
#'@param na.rm logical: should missing values (including \code{NA} and \code{NaN})
#'be omitted from the calculations? Default is \code{FALSE}.
#'
#'@return A list with the following elements when the set p-value is computed :\itemize{
#'   \item \code{set_score_obs}: the approximation of the observed set score
#'   \item \code{set_pval}: the associated set p-value
#' }
#' or a list with the following elements when gene-wise p-values are computed:\itemize{
#'   \item \code{gene_scores_obs}: vector of approximating the observed gene-wise scores
#'   \item \code{gene_pvals}: vector of associated gene-wise p-values
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
#'t <- matrix(rep(1:(r/4)), 4, ncol=1, nrow=r)
#'sigma <- 0.4
#'b0 <- 1
#'
#'#under the null:
#'b1 <- 0
#'#under the alternative:
#'#b1 <- 0.5
#'y.tilde <- b0 + b1*t + rnorm(r, sd = sigma)
#'y <- t(matrix(rnorm(n*r, sd = sqrt(sigma*abs(y.tilde))), ncol=n, nrow=r) +
#'       matrix(rep(y.tilde, n), ncol=n, nrow=r))
#'x <- matrix(1, ncol=1, nrow=r)
#'
#'#run test
#'asymTestRes <- vc_test_asym(y, x, phi=cbind(t, t^2), w=matrix(1, ncol=ncol(y), nrow=nrow(y)),
#'                            Sigma_xi=diag(2), indiv=1:r, genewise_pvals=TRUE)
#'plot(density(asymTestRes$gene_pvals))
#'quantile(asymTestRes$gene_pvals)
#'
#'@importFrom CompQuadForm davies
#'@importFrom stats pchisq var
#'
#'@export
vc_test_asym <- function(y, x, indiv=rep(1,nrow(x)), phi, w, Sigma_xi = diag(ncol(phi)),
                         genewise_pvals=FALSE, homogen_traj=FALSE, na.rm = FALSE){

  if(homogen_traj){
    score_list <- vc_score_h(y = y, x = x, indiv = factor(indiv), phi = phi, w = w,
                             Sigma_xi = Sigma_xi, na_rm = na.rm)
  }else{
    score_list <- vc_score(y = y, x = x, indiv = factor(indiv), phi = phi, w = w,
                           Sigma_xi = Sigma_xi, na_rm = na.rm)
  }

  nindiv <- nrow(score_list$q_ext)
  ng <- nrow(y)
  nphi <- ncol(phi)

  if(ng*nindiv < 1){
    stop("no gene measured/no sample included ...")
  }

  if(genewise_pvals){
    gene_scores_obs <- score_list$gene_scores_unscaled
    if(nindiv == 1){
      pv <- stats::pchisq(gene_scores_obs, df = 1, lower.tail = FALSE)
    }else if(nphi == 1){
      gene_lambda <- apply(X=score_list$q_ext, MARGIN = 2, FUN = stats::var)
      if(ng == 1){
        pv <- stats::pchisq(gene_scores_obs/gene_lambda, df = 1, lower.tail = FALSE)
      }else{
        pv <- unlist(mapply(FUN=CompQuadForm::davies, q=gene_scores_obs, lambda=gene_lambda, lim=15000, acc=0.0005)["Qq",])
        #pv <- stats::pchisq(gene_scores_obs/gene_lambda, df = 1, lower.tail = FALSE) # same result ? only if phi is univariate
      }
    }else{
      gene_inds <- lapply(1:ng, function(x){x + (ng)*(0:(nphi-1))})

      gene_lambda <- lapply(gene_inds, function(x){
        Sig_q_gene <- cov(score_list$q_ext[, x, drop=FALSE])
        lam <- try(svd(Sig_q_gene)$d)
        if (inherits(lam, "try-error")){
          lam <- try(svd(round(Sig_q_gene, 6))$d)
          if (inherits(lam, "try-error")){
            warning("Error in svd decomposition for at least one gene")
            lam <- NA
          }
        }
        return(lam)
      })

      pv <- unlist(mapply(FUN=CompQuadForm::davies, q=gene_scores_obs, lambda=gene_lambda, lim=15000, acc=0.0005)["Qq",])
    }

    names(pv) <- rownames(y)
    ans <- list("gene_scores_obs" = gene_scores_obs, "gene_pvals" = pv)

  }else{

    if(nindiv == 1){
      Sig_q <- matrix(1, ng, ng)
    }else{
      Sig_q <- cov(score_list$q_ext)
    }

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

