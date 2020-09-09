#'Computes variance component test statistic and its permuted distribution
#'
#'This function computes an approximation of the Variance Component test for a
#'mixture of \eqn{\chi^{2}}s using Davies method from \code{\link[CompQuadForm]{davies}}
#'
#'@keywords internal
#'
#'@param y a numeric matrix of dim \code{g x n} containing the raw RNA-seq counts for g
#'genes from \code{n} samples
#'
#'@param x a numeric design matrix of dim \code{n x p} containing the \code{p} covariates
#' to be adjusted for
#'
#'@param indiv a vector of length \code{n} containing the information for
#'attributing each sample to one of the studied individuals. Coerced
#'to be a \code{factor}
#'
#'@param phi a numeric design matrix of size \code{n x K} containing the \code{K} variables
#'to be tested
#'
#'@param w a vector of length \code{n} containing the weights for the \code{n}
#'samples.
#'
#'@param Sigma_xi a matrix of size \code{K x K} containing the covariance matrix
#'of the \code{K} random effects on \code{phi}
#'
#'@param na_rm logical: should missing values (including \code{NA} and \code{NaN})
#'be omitted from the calculations? Default is \code{FALSE}.
#'
#'@param n_perm the number of permutation to perform. Default is \code{1000}.
#'
#'@param progressbar logical indicating whether a progress bar should be displayed
#'when computing permutations (only in interactive mode).
#'
#'@param parallel_comp a logical flag indicating whether parallel computation
#'should be enabled. Only Linux and MacOS are supported, this is ignored on Windows.
#'Default is \code{TRUE}.
#'
#'@param nb_cores an integer indicating the number of cores to be used when
#'\code{parallel_comp} is \code{TRUE}.
#'Only Linux and MacOS are supported, this is ignored on Windows.
#'Default is \code{parallel::detectCores() - 1}.
#'
#'@return A list with the following elements:\itemize{
#'   \item \code{score}: an approximation of the observed set score
#'   \item \code{scores_perm}: a vector containing the permuted set scores
#'   \item \code{gene_scores_unscaled}: approximation of the individual gene scores
#'   \item \code{gene_scores_unscaled_perm}: a list of approximation of the permuted individual gene scores
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
#'t <- matrix(rep(1:3), r/3, ncol=1, nrow=r)
#'sigma <- 0.4
#'b0 <- 1
#'
#'#under the null:
#'b1 <- 0
#'#under the alternative:
#'b1 <- 0.7
#'y.tilde <- b0 + b1*t + rnorm(r, sd = sigma)
#'y <- t(matrix(rnorm(n*r, sd = sqrt(sigma*abs(y.tilde))), ncol=n, nrow=r) +
#'       matrix(rep(y.tilde, n), ncol=n, nrow=r))
#'x <- matrix(1, ncol=1, nrow=r)
#'
#'#run test
#'scoreTest <- vc_score_perm(y, x, phi=t, w=matrix(1, ncol=ncol(y), nrow=nrow(y)),
#'                     Sigma_xi=matrix(1), indiv=rep(1:(r/3), each=3), parallel_comp = FALSE)
#'scoreTest$score
#'
#'@importFrom CompQuadForm davies
#'@importFrom stats model.matrix
#'@importFrom pbapply pbsapply
#'@importFrom parallel mclapply
#'
#'@export
vc_score_perm <- function(y, x, indiv, phi, w, Sigma_xi = diag(ncol(phi)),
                          na_rm = FALSE, n_perm = 1000, progressbar = TRUE,
                          parallel_comp = TRUE, nb_cores = parallel::detectCores() - 1) {
  ## validity checks
  if(sum(!is.finite(w))>0){
    stop("At least 1 non-finite weight in 'w'")
  }

  ## dimensions check------
  stopifnot(is.matrix(y))
  stopifnot(is.matrix(x))
  stopifnot(is.matrix(phi))

  g <- nrow(y) # the number of genes measured
  n <- ncol(y) # the number of samples measured
  p <- ncol(x) # the number of covariates
  n_t <- ncol(phi) # the number of time bases
  stopifnot(nrow(x) == n)
  stopifnot(nrow(w) == g)
  stopifnot(ncol(w) == n)
  stopifnot(nrow(phi) == n)
  stopifnot(length(indiv) == n)


  # the number of random effects
  if (length(Sigma_xi) == 1) {
    K <- 1
    Sigma_xi <- matrix(Sigma_xi, K, K)
  } else {
    K <- nrow(Sigma_xi)
    stopifnot(ncol(Sigma_xi) == K)
  }
  stopifnot(n_t == K)


  ## data formating ------
  indiv <- as.factor(indiv)
  nb_indiv <- length(levels(indiv))

  ## OLS for conditional mean -----
  y_T <- t(y)
  if(na_rm & sum(is.na(y_T))>0){
    y_T0 <- y_T
    y_T0[is.na(y_T0)] <- 0
    yt_mu <- y_T - x%*%solve(crossprod(x))%*%t(x)%*%y_T0
    rm(y_T0)
  }else{
    yt_mu <- y_T - x%*%solve(crossprod(x))%*%t(x)%*%y_T
  }
  rm(y_T, y)
  # x_tilde_list <- y_tilde_list <- Phi_list <- list()
  # for (i in 1:nb_indiv) {
  #   select <- indiv==levels(indiv)[i]
  #   n_i <- length(which(select))
  #   x_i <- x[select,]
  #   y_i <- y[,select]
  #   phi_i <- phi[select,]
  #   Phi_list[[i]] <- kronecker(diag(g), phi_i)
  #   x_tilde_list[[i]] <- kronecker(diag(g), x_i)
  #   y_tilde_list[[i]] <- matrix(t(y_i), ncol=1)
  # }
  # x_tilde <- do.call(rbind, x_tilde_list)
  # y_tilde <- do.call(rbind, y_tilde_list)
  # Phi <- do.call(rbind, Phi_list)
  #
  # alpha <- solve(t(x_tilde)%*%x_tilde)%*%t(x_tilde)%*%y_tilde
  # mu_new <- x_tilde %*% alpha
  # y_mu <- y_tilde - mu_new


  ## test statistic computation ------
  # q <- matrix(NA, nrow=nb_indiv, ncol=g*K)
  # XT_i <- array(NA, c(nb_indiv, g*p, g*K))
  # U <- matrix(NA, nrow = nb_indiv, ncol = p*g)
  #
  # long_indiv <- rep(indiv, each = g)
  # xtx_inv <- solve(t(x_tilde) %*% x_tilde)
  # Sigma_xi_new_sqrt <- kronecker(diag(g), (Sigma_xi*diag(K))%^% (-0.5))
  #
  # for (i in 1:nb_indiv){
  #   #for all the genes at once
  #   select <- indiv==levels(indiv)[i]
  #   long_select <- long_indiv==levels(indiv)[i]
  #   y_mu_i <-  as.vector(y_mu[long_select,])
  #   # y_tilde_i <- c(t(y_ij))
  #   x_tilde_i <- x_tilde[long_select,]
  #
  #   sigma_eps_inv_diag <- as.vector(t(w)[select,])#/sigma
  #   T_i <- sigma_eps_inv_diag*(Phi[long_select,] %*% Sigma_xi_new_sqrt)
  #   q[i,] <- c(y_mu_i %*% T_i)
  #   XT_i[i,,] <- t(x_tilde_i) %*% T_i
  #   U[i,] <- xtx_inv %*% t(x_tilde_i) %*% y_mu_i
  # }
  # XT <- colMeans(XT_i)
  # q_ext <- q - U %*% XT

  #sig_eps_inv <- w/sigma #no need as this just scales the test statistics
  # browser()
  sig_xi_sqrt <- (Sigma_xi*diag(K))%^% (-0.5)
  sig_eps_inv_T <- t(w)

  if(length(levels(indiv))>1){
    indiv_mat <- stats::model.matrix(~0 + factor(indiv))
  }else{
    indiv_mat <- matrix(as.numeric(indiv), ncol=1)
  }
  avg_xtx_inv_tx <- nb_indiv*tcrossprod(solve(crossprod(x, x)), x)

  compute_genewise_scores <- function(v, indiv_mat, avg_xtx_inv_tx){
    phi_perm <- phi[v, , drop=FALSE]
    phi_sig_xi_sqrt <- phi_perm%*%sig_xi_sqrt
    T_fast <- do.call(cbind, replicate(K, sig_eps_inv_T, simplify = FALSE))*matrix(apply(phi_sig_xi_sqrt, 2, rep, g), ncol = g*K)
    q_fast <- matrix(yt_mu, ncol=g*n_t, nrow=n)*T_fast
    if(na_rm & sum(is.na(q_fast))>0){
      q_fast[is.na(q_fast)] <- 0
    }
    q <- crossprod(indiv_mat, q_fast)
    XT_fast <- t(x)%*%T_fast/nb_indiv
    U_XT <- matrix(yt_mu, ncol=g*n_t, nrow=n)*crossprod(avg_xtx_inv_tx, XT_fast)
    if(na_rm & sum(is.na(U_XT))>0){
      U_XT[is.na(U_XT)] <- 0
    }
    U_XT_indiv <- crossprod(indiv_mat, U_XT)
    q_ext <-  q - U_XT_indiv
    qq <- colSums(q, na.rm = na_rm)^2/nb_indiv
    return(rowSums(matrix(qq, ncol=K)))# genewise scores
  }

  o <- order(as.numeric(unlist(split(x = as.character(1:n), f = indiv))))
  perm_list <- c(list(1:n), lapply(1:n_perm, function(x){as.numeric(unlist(lapply(split(x = as.character(1:n), f = indiv), FUN=sample)))[o]}))

  if(!parallel_comp){
    if(progressbar){
      gene_Q <- pbapply::pbsapply(X = perm_list, FUN = compute_genewise_scores,
                                  indiv_mat = indiv_mat, avg_xtx_inv_tx = avg_xtx_inv_tx)
    }else{
      gene_Q <- sapply(X = perm_list, FUN = compute_genewise_scores,
                       indiv_mat = indiv_mat, avg_xtx_inv_tx = avg_xtx_inv_tx)
    }
  }else{
    if(progressbar){
      gene_Q <- pbapply::pbsapply(perm_list, compute_genewise_scores,
                                  indiv_mat = indiv_mat, avg_xtx_inv_tx = avg_xtx_inv_tx,
                                  cl = nb_cores)
    }else{
      gene_Q <- simplify2array(parallel::mclapply(X = perm_list, FUN = compute_genewise_scores,
                                                  indiv_mat = indiv_mat, avg_xtx_inv_tx = avg_xtx_inv_tx,
                                                  mc.cores = nb_cores))
    }
  }
  rownames(gene_Q) <- colnames(yt_mu)
  QQ <- colSums(gene_Q)

  return(list("score" = QQ[1], "scores_perm" = QQ[-1],
              "gene_scores_unscaled" = gene_Q[,1], "gene_scores_unscaled_perm" = gene_Q[,-1]))
}
