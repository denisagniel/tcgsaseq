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
#'ng <- 100
#'nindiv <- 30
#'nt <- 5
#'nsample <- nindiv*nt
#'tim <- matrix(rep(1:nt), nindiv, ncol=1, nrow=nsample)
#'tim <- cbind(tim, tim^2)
#'sigma <- 5
#'b0 <- 10
#'
#'#under the null:
#'beta1 <- rnorm(n=ng, 0, sd=0)
#'#under the (heterogen) alternative:
#'beta1 <- rnorm(n=ng, 0, sd=0.1)
#'#under the (homogen) alternative:
#'beta1 <- rnorm(n=ng, 0.06, sd=0)
#'
#'y.tilde <- b0 + rnorm(ng, sd = sigma)
#'y <- t(matrix(rep(y.tilde, nsample), ncol=ng, nrow=nsample, byrow=TRUE) +
#'       matrix(rep(beta1, each=nsample), ncol=ng, nrow=nsample, byrow=FALSE)*matrix(rep(tim, ng),
#'                                                             ncol=ng, nrow=nsample, byrow=FALSE) +
#'       matrix(rnorm(ng*nsample, sd = sigma), ncol=ng, nrow=nsample, byrow=FALSE)
#'       )
#'myindiv <- rep(1:nindiv, each=nt)
#'x <- cbind(1, myindiv/2==floor(myindiv/2))
#'myw <- matrix(rnorm(nsample*ng, sd=0.1), ncol=nsample, nrow=ng)
#'
#'#run test
#'score_homogen <- vc_score_h(y, x, phi=tim, indiv=myindiv,
#'                            w=myw, Sigma_xi=cov(tim))
#'
#'score_heterogen <- vc_score(y, x, phi=tim, indiv=myindiv,
#'                            w=myw, Sigma_xi=cov(tim))
#'
#'scoreTest_homogen <- vc_test_asym(y, x, phi=tim, indiv=rep(1:nindiv, each=nt),
#'                                  w=matrix(1, ncol=ncol(y), nrow=nrow(y)), Sigma_xi=cov(tim),
#'                                  homogen_traj = TRUE)
#'scoreTest_homogen$set_pval
#'scoreTest_heterogen <- vc_test_asym(y, x, phi=tim, indiv=rep(1:nindiv, each=nt),
#'                                    w=matrix(1, ncol=ncol(y), nrow=nrow(y)), Sigma_xi=cov(tim),
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

  # x_tilde_list <- y_tilde_list <- Phi_list <- list()
  # for (i in 1:nb_indiv) {
  #   select <- indiv==levels(indiv)[i]
  #   n_i <- length(which(select))
  #   x_i <- x[select,]
  #   y_i <- y_T[select,]
  #   phi_i <- phi[select,]
  #   Phi_list[[i]] <- do.call(rbind, replicate(g, phi_i, simplify = FALSE)) #TODO
  #   x_tilde_list[[i]] <- matrix(data=rep(x_i, each=g), ncol = p) #TODO
  #   y_tilde_list[[i]] <- matrix(y_i, ncol=1)
  # }
  # x_tilde <- do.call(rbind, x_tilde_list)
  # y_tilde <- do.call(rbind, y_tilde_list)
  # Phi <- do.call(rbind, Phi_list)
  #
  # alpha <- solve(t(x_tilde)%*%x_tilde)%*%t(x_tilde)%*%y_tilde
  # mu_new <- x_tilde %*% alpha
  # y_mu <- y_tilde - mu_new


  alpha <- solve(crossprod(x))%*%t(x)%*%rowMeans(y_T)
  yt_mu <- y_T - do.call(cbind,replicate(g, x%*%alpha, simplify = FALSE))


  ## test statistic computation ------
  sig_xi_sqrt <- (Sigma_xi*diag(K))%^% (-0.5) #sig_xi_sqrt <- (Sigma_xi %^% (-0.5))
  # xtx_inv <- solve(t(x_tilde) %*% x_tilde)
  # long_indiv <- rep(indiv, each = g)
  #
  # q <- matrix(NA, nrow=nb_indiv, ncol=K)
  # XT_i <- array(NA, c(nb_indiv, p, K))
  # U <- matrix(NA, nrow = nb_indiv, ncol = p)
  #
  # for (i in 1:nb_indiv){
  #   #for all the genes at once
  #   select <- indiv==levels(indiv)[i]
  #   long_select <- long_indiv==levels(indiv)[i]
  #   y_mu_i <-  as.vector(y_mu[long_select,])
  #   # y_tilde_i <- c(y_ij)
  #   x_tilde_i <- x_tilde[long_select,]
  #
  #   sigma_eps_inv_diag <- c(t(w)[select,])
  #   T_i <- sigma_eps_inv_diag*(Phi[long_select,] %*% sig_xi_sqrt)
  #   q[i,] <- c(y_mu_i %*% T_i)
  #   XT_i[i,,] <- t(x_tilde_i) %*% T_i
  #   U[i,] <- xtx_inv %*% t(x_tilde_i) %*% y_mu_i
  # }
  # XT <- colMeans(XT_i)
  # q_ext <- q - U %*% XT

  sig_eps_inv_T <- t(w)
  phi_sig_xi_sqrt <- phi%*%sig_xi_sqrt
  T_fast <- do.call(cbind, replicate(K, sig_eps_inv_T, simplify = FALSE))*t(matrix(rep(t(phi_sig_xi_sqrt), each=g), nrow=g*K))
  q_fast <- do.call(cbind, replicate(K, yt_mu, simplify = FALSE))*T_fast
  q <- do.call(rbind, by(do.call(cbind, by(t(q_fast), as.factor(rep(1:K, each=g)), FUN=colSums, simplify=FALSE)), indiv, FUN=colSums, simplify=FALSE))

  avg_xtx_inv_tx <- solve(t(x)%*%x)%*%t(x)
  #XT <- matrix(rowSums(crossprod(x, T_fast))/nb_indiv, nrow=p)
  XT_fast <- crossprod(x, T_fast)/nb_indiv
  U_XT_fast <- rowMeans(yt_mu)*rowSums(t(avg_xtx_inv_tx)%*%XT_fast)
  U_XT_fast <- rowMeans(yt_mu)*do.call(cbind, by(crossprod(XT_fast, avg_xtx_inv_tx), rep(1:K, each=g), FUN=colSums, simplify=FALSE))
  U_XT_indiv <- do.call(rbind, by(U_XT_fast, indiv, FUN=colSums, simplify=FALSE))
  #U2 <- apply(avg_xtx_inv_tx, 1, function(xx){rowMeans(do.call(cbind, replicate(g, xx, simplify = FALSE))*yt_mu)})
  #U_indiv <- do.call(rbind, by(U2, indiv, FUN=colSums, simplify=FALSE))
  q_ext <-  q - U_XT_indiv #U_indiv %*% XT

  qq <- colSums(q)^2/nb_indiv
  #do.call(rbind, by(t(qq), rep(1:K, g), FUN=colSums, simplify=FALSE))
  QQ <- sum(qq)

  return(list("score"=QQ, "q" = q, "q_ext"=q_ext,
              "gene_scores_unscaled"=qq))
}
