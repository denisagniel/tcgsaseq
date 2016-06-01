#'Computes variance component test statistic for homogeneous trajectory
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
#'   \item \code{score}: approximation of the observed score
#'   \item \code{q}: TODO
#'   \item \code{q_ext}: TODO
#' }
#'
#'@seealso \code{\link[CompQuadForm]{davies}}
#'
#'@examples
#'
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
  qq <- ncol(x) # the number of covariates
  n.t <- ncol(phi) # the number of time bases
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
  stopifnot(n.t == K)
  
  
  ## data formating ------
  indiv <- as.factor(indiv)
  nb_indiv <- length(levels(indiv))
  
  x_tilde_list <- y_tilde_list <- Phi_list <- list()
  for (i in 1:nb_indiv) {
    select <- indiv==levels(indiv)[i]
    n_i <- length(which(select))
    x_i <- x[select,]
    y_i <- y.t[,select]
    phi_i <- phi[select,]
    Phi_list[[i]] <- rep(phi_i,g) %>% matrix(ncol = n.t)
    x_tilde_list[[i]] <- rep(x_i,each=g) %>% matrix(ncol = qq)
    y_tilde_list[[i]] <- matrix(t(y_i), ncol=1)
  }
  x_tilde <- do.call(rbind, x_tilde_list)
  y_tilde <- do.call(rbind, y_tilde_list)
  Phi <- do.call(rbind, Phi_list)
  
  alpha <- solve(t(x_tilde)%*%x_tilde)%*%t(x_tilde)%*%y_tilde
  mu_new <- x_tilde %*% alpha
  y_mu <- y_tilde - mu_new
  
  xtx_inv <- solve(t(x_tilde) %*% x_tilde)
  Sigma_xi_sqrt <- (Sigma_xi %^% (-0.5))
  #browser()
  
  ## test statistic computation ------
  #q <- matrix(NA, nrow=nb_indiv, ncol=K)
  q <- matrix(NA, nrow=nb_indiv, ncol=K)
  XT_i <- array(NA, c(nb_indiv, qq, K))
  U <- matrix(NA, nrow = nb_indiv, ncol = qq)
  
  long_indiv <- rep(indiv, each = g)
  
  for (i in 1:nb_indiv){
    select <- indiv==levels(indiv)[i]
    long_select <- long_indiv==levels(indiv)[i]
    y_mu_i <-  as.vector(y_mu[select,])
    x_tilde_i <- x_tilde[long_select,]
    
    sigma_eps_inv_diag <- c(t(w.t[,select]))
    T_i <- sigma_eps_inv_diag*(Phi[long_select,] %*% Sigma_xi_sqrt)
    q[i,] <- c(y_mu_i %*% T_i)
    XT_i[i,,] <- t(x_tilde_i) %*% T_i
    U[i,] <- xtx_inv %*% t(x_tilde_i) %*% y_mu_i
  }
  XT <- colMeans(XT_i)
  q_ext <- q - U %*% XT
  QQ <- sum(colSums(q)^2/nrow(q))
  
  return(list("score"=QQ, "q" = q, "q_ext"=q_ext))
}