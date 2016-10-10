#'Computes variance component test statistic for homogeneous trajectory
#'
#'This function computes an approximation of the variance component test for homogeneous trajectory
#'based on the asymptotic distribution of a mixture of \eqn{\chi^{2}}s using Davies method
#'from \code{\link[CompQuadForm]{davies}}
#'
#'@keywords internal
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
#'@return A list with the following elements:\itemize{
#'   \item \code{score}: approximation of the set observed score
#'   \item \code{q}: observation-level contributions to the score
#'   \item \code{q_ext}: psuedo-observations used to compute covariance taking into account the contributions of OLS estimates
#'   \item \code{gene_scores}: approximation of the individual gene scores
#' }
#'
#'
#'@examples
#'rm(list=ls())
#'set.seed(123)
#'
#'##generate some fake data
#'########################
#'n <- 100
#'nindiv <- 100
#'nt <- 30
#'r <- nindiv*nt
#'t <- matrix(rep(1:nt), nindiv, ncol=1, nrow=r)
#'sigma <- 0.4
#'b0 <- 1
#'
#'#under the null:
#'beta1 <- rep(rnorm(n=nindiv, 0, sd=0), each=nt)
#'#under the alternative:
#'beta1 <- rep(rnorm(n=nindiv, 0, sd=0.1), each=nt)
#'
#'y.tilde <- b0 + beta1*t + rnorm(r, sd = sigma)
#'y <- t(matrix(rnorm(n*r, sd = sqrt(sigma*abs(y.tilde))), ncol=n, nrow=r) +
#'       matrix(rep(y.tilde, n), ncol=n, nrow=r))
#'x <- matrix(1, ncol=1, nrow=r)
#'
#'#run test
#'score_homogen <- vc_score_h(y, x, phi=t, indiv=rep(1:nindiv, each=nt),
#'                            w=matrix(1, ncol=ncol(y), nrow=nrow(y)),
#'                            Sigma_xi=matrix(1))
#'score_heterogen <- vc_score(y, x, phi=t, indiv=rep(1:nindiv, each=nt),
#'                            w=matrix(1, ncol=ncol(y), nrow=nrow(y)),
#'                            Sigma_xi=matrix(1))
#'scoreTest_homogen <- vc_test_asym(y, x, phi=t, indiv=rep(1:nindiv, each=nt),
#'                                  w=matrix(1, ncol=ncol(y), nrow=nrow(y)), Sigma_xi=matrix(1),
#'                                  homogen_traj = TRUE)
#'scoreTest_homogen$set_pval
#'scoreTest_heterogen <- vc_test_asym(y, x, phi=t, indiv=rep(1:nindiv, each=nt),
#'                                    w=matrix(1, ncol=ncol(y), nrow=nrow(y)), Sigma_xi=matrix(1),
#'                                    homogen_traj = FALSE)
#'scoreTest_heterogen$set_pval
#'
#'@seealso \code{\link[CompQuadForm]{davies}}
#'@importFrom CompQuadForm davies
#'
#'@export
vc_score_h <- function(y, x, indiv, phi, w, Sigma_xi = diag(ncol(phi))) {

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

  y_T <- t(y)

  x_tilde_list <- y_tilde_list <- Phi_list <- list()
  for (i in 1:nb_indiv) {
    select <- indiv==levels(indiv)[i]
    n_i <- length(which(select))
    x_i <- x[select,]
    y_i <- y_T[select,]
    phi_i <- phi[select,]
    Phi_list[[i]] <- matrix(rep(phi_i, g), ncol = n_t) #TODO
    x_tilde_list[[i]] <- matrix(data=rep(x_i, each=g), ncol = p) #TODO
    y_tilde_list[[i]] <- matrix(y_i, ncol=1)
  }
  x_tilde <- do.call(rbind, x_tilde_list)
  y_tilde <- do.call(rbind, y_tilde_list)
  Phi <- do.call(rbind, Phi_list)

  alpha <- solve(t(x_tilde)%*%x_tilde)%*%t(x_tilde)%*%y_tilde
  mu_new <- x_tilde %*% alpha
  y_mu <- y_tilde - mu_new

  xtx_inv <- solve(t(x_tilde) %*% x_tilde)
  Sigma_xi_sqrt <- (Sigma_xi %^% (-0.5)) #TODO

  ## test statistic computation ------
  q <- matrix(NA, nrow=nb_indiv, ncol=K)
  XT_i <- array(NA, c(nb_indiv, p, K))
  U <- matrix(NA, nrow = nb_indiv, ncol = p)

  long_indiv <- rep(indiv, each = g)

  for (i in 1:nb_indiv){
    #for all the genes at once
    select <- indiv==levels(indiv)[i]
    long_select <- long_indiv==levels(indiv)[i]
    y_mu_i <-  as.vector(y_mu[long_select,])
    # y_tilde_i <- c(y_ij)
    x_tilde_i <- x_tilde[long_select,]

    sigma_eps_inv_diag <- c(w[,select])
    T_i <- sigma_eps_inv_diag*(Phi[long_select,] %*% Sigma_xi_sqrt)
    q[i,] <- c(y_mu_i %*% T_i)
    XT_i[i,,] <- t(x_tilde_i) %*% T_i
    U[i,] <- xtx_inv %*% t(x_tilde_i) %*% y_mu_i
  }
  XT <- colMeans(XT_i)
  q_ext <- q - U %*% XT
  qq <- colSums(q)^2/nb_indiv
  QQ <- sum(qq)

  return(list("score"=QQ, "q" = q, "q_ext"=q_ext,
              "gene_scores_unscaled"=qq))
}
