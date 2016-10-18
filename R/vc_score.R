#'Computes variance component test statistic for longitudinal
#'
#'This function computes an approximation of the Variance Component test for a
#'mixture of \eqn{\chi^{2}}s using Davies method from \code{\link[CompQuadForm]{davies}}
#'
#'@keywords internal
#'
#'@param y a numeric matrix of dim \code{g x n} containing the raw RNAseq counts for g
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
#'@return A list with the following elements:\itemize{
#'   \item \code{score}: approximation of the set observed score
#'   \item \code{q}: observation-level contributions to the score
#'   \item \code{q_ext}: psuedo-observations used to compute covariance taking into account the contributions of OLS estimates
#'   \item \code{gene_scores}: approximation of the individual gene scores
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
#'b1 <- 0.7
#'y.tilde <- b0 + b1*t + rnorm(r, sd = sigma)
#'y <- t(matrix(rnorm(n*r, sd = sqrt(sigma*abs(y.tilde))), ncol=n, nrow=r) +
#'       matrix(rep(y.tilde, n), ncol=n, nrow=r))
#'x <- matrix(1, ncol=1, nrow=r)
#'
#'#run test
#'scoreTest <- vc_score(y, x, phi=t, w=matrix(1, ncol=ncol(y), nrow=nrow(y)),
#'                     Sigma_xi=matrix(1), indiv=rep(1:4, each=3))
#'scoreTest$score
#'
#'@importFrom CompQuadForm davies
#'
#'@export
vc_score <- function(y, x, indiv, phi, w, Sigma_xi = diag(ncol(phi))) {
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
  yt_mu <- y_T - x%*%solve(crossprod(x))%*%t(x)%*%y_T

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
  sig_xi_sqrt <- (Sigma_xi*diag(K))%^% (-0.5)
  sig_eps_inv_T <- t(w)
  phi_sig_xi_sqrt <- phi%*%sig_xi_sqrt
  T_fast <- do.call(cbind, replicate(K, sig_eps_inv_T, simplify = FALSE))*t(matrix(rep(t(phi_sig_xi_sqrt), each=g), nrow=g*K))
  q_fast <- matrix(yt_mu, ncol=g*n_t, nrow=n)*T_fast

  q <- do.call(rbind, by(q_fast, indiv, FUN=colSums, simplify=FALSE))
  XT_fast <- t(x)%*%T_fast/nb_indiv
  avg_xtx_inv_tx <- nb_indiv*solve(t(x)%*%x)%*%t(x)
  U_XT <- matrix(yt_mu, ncol=g*n_t, nrow=n)*t(avg_xtx_inv_tx)%*%XT_fast
  U_XT_indiv <- do.call(rbind, by(U_XT, indiv, FUN=colSums, simplify=FALSE))
  q_ext <-  q - U_XT_indiv
  #sapply(1:6, function(i){(q_ext[i,] - q_ext_fast_indiv[i,])})

  qq <- colSums(q)^2/nb_indiv # genewise scores ?
  #do.call(rbind, by(t(qq), rep(1:K, g), FUN=colSums, simplify=FALSE))
  QQ <- sum(qq) #nb_indiv=nrow(q) # set score

  return(list("score"=QQ, "q" = q, "q_ext"=q_ext,
              "gene_scores_unscaled" = qq)) #gene_score(qq, q_ext)
}

gene_score <- function(indiv_q_2, q_ext) {
  indiv_v <- apply(X=q_ext, margin = 2, FUN=stats::var)
  indiv_chi <- indiv_q_2/indiv_v
}
