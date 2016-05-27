
#' Computing negative binomial simulations results
#'
#'@examples
#'
#'\dontrun{
#'
#'all_sim_pvals <- replicate(1000, nb_sim_fn(type = 'standard_null'))
#'}
#'
#'@keywords internal
#'@importFrom stats rnorm qnorm pnorm integrate
#@importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#@importFrom edgeR DGEList estimateDisp calcNormFactors
#@importFrom limma roast
#'@export


nb_sim_fn <- function(n = 250,
                      n_t = 5,
                      beta = 2,
                      nGenes = 1000,
                      type) {
  sim_data <- sim_nb_data(n = n, n_t = n_t,
                          beta = beta, nGenes = nGenes, type = type)

  set_size <- 10
  set_ind <- 1:set_size
  groups <- nGenes/set_size %>% floor
  
  x <- sim_data$x
  y <- sim_data$y
  tt <- sim_data$tt
  indiv <- sim_data$indiv
  
  design_r <- cbind(x, tt)
  
  w <- sp_weights(x = x, 
                            y = y, 
                            phi = matrix(tt, ncol = 1), 
                            preprocessed = TRUE, 
                            doPlot = FALSE)
  voom_w <- voom_weights(x = design_r, 
                                   y = t(y), 
                                   preprocessed = TRUE, 
                                   doPlot = FALSE)
  
  genesets <- data.frame(Symbol = rep(1:groups, each = set_size),
                         Chr = rep(1:groups, each = set_size))
  
  ydge <- edgeR::DGEList(counts=t(y), genes=genesets)
  ydge <- edgeR::calcNormFactors(ydge)
  ydge <- edgeR::estimateDisp(ydge, design_r, robust=TRUE)
  
  y_set <- y[,set_ind]
  w_set <- w[set_ind,]
  voomw_set <- voom_w[set_ind,]
  
  
  index <- list(X=(1:nGenes %in% set_ind))
  tcg <- vc_test_asym(y = t(y_set), 
                                x = x, 
                                indiv = indiv, 
                                phi = matrix(tt), 
                                Sigma_xi = 1, 
                                w = w_set)
  
  roast_voom <- limma::roast(y=t(y_set), 
                             set.statistic = 'msq',
                             design=design_r, 
                             contrast=3, 
                             block = indiv, 
                             correlation = 0,
                             weights = voomw_set)
  roast_edger <- limma::roast(y=ydge, 
                              index = index, 
                              set.statistic = 'msq',
                              design=design.r, 
                              contrast=3, 
                              block = indiv, 
                              correlation = 0)
  deseq <- deseq_fn(y_set)
  pvals <- data.frame(
    'tcgsaseq'=tcg$pval,
    'roast_voom'=roast_voom$p.value["Mixed", "P.Value"],
    'roast_edger'=roast_edger$PValue.Mixed,
    'deseq'=deseq
  )
  pvals
}

sim_nb_data <- function(n = 250,
                        n_t = 5,
                        beta = 2,
                        nGenes = 1000,
                        type) {
  stopifnot(type %in% c('standard_null', 'alternative', 'td_variance')) 
  mu_x <- rexp(n)*10
  x <- cbind(1, rep(rnorm(n, mean = mu_x), each = n_t))
  tt <- c(replicate(n, sample(1:10, size = n_t)))
  b_0 <- replicate(nGenes, rep(rnorm(n), each = n_t))
  if (type == 'standard_null') {
    y <- replicate(nGenes, rnbinom(n*n_t, mu = 100 + b_0 + rowSums(x), size = 1)) + 1
  } else if (type == 'td_variance') {
    y <- replicate(nGenes, rnbinom(n*n_t, mu = 100 + b_0 + rowSums(x), size = 1/tt)) + 1
  } else if (type == 'alternative') {
    y <- cbind(replicate(2, rnbinom(n*n_t, mu = tt*beta + 100 + b_0 + rowSums(x), size = 1)),
               replicate(nGenes-2, rnbinom(n*n_t, mu = 100 + b_0 + rowSums(x), size = 1))) + 1
  }
  
  indiv <- rep(1:n, each = n_t)
  
  list(x = x, tt = tt, y = y, indiv = indiv)
}


Kappa_j <- function(r, alpha_DEseq){
  sig <- stats::qnorm(1-alpha_DEseq/2)
  temp_int <- stats::integrate(function(x){exp(-x^2/2)*stats::pnorm((r*x-sig)/sqrt(1-r^2))},
                               lower=-sig,
                               upper=sig)
  1/log(1-alpha_DEseq)*log(1-(1/(1-alpha_DEseq)*sqrt(2/pi)*temp_int$value))
}

deseq_fn <- function(y) {
  # browser()
  y_dsq <- DESeq2::DESeqDataSetFromMatrix(countData = t(y),
                                          colData = cbind.data.frame("indiv"=as.factor(indiv),
                                                                     "time"=as.numeric(tt),
                                                                     "x"=as.numeric(x[,2])),
                                          design = ~ x + time)
  res_dsq <- DESeq2::DESeq(y_dsq, test="LRT", reduced = ~ x)
  pvals <- DESeq2::results(res_dsq)$pvalue
  
  pmin <- min(pvals, na.rm = TRUE)
  Rij <- cor(y)
  Rj <- sapply(1:ncol(Rij), function(j){max(abs(Rij[1:(j-1),j]), na.rm = TRUE)})
  if(length(which(is.infinite(Rj)))>0){
    Rj[which(is.infinite(Rj))] <- 0
  }
  
  Keff <- 1 + sum(sapply(Rj[-1], Kappa_j, alpha_DEseq = 0.05))
  Keff_approx <- 1 + sum(sqrt(1-Rj[-1]^(-1.31*log10(0.05))))
  Meff <- 1 + sum(1-cor(y)^2, na.rm = TRUE)/ncol(Rij) #Cheverud–Nyholt
  cbind("minTest_exact" = 1-(1-pmin)^Keff,
        "minTest_approx" = 1-(1-pmin)^Keff_approx,
        "minTest_CN" = 1-(1-pmin)^Meff)
  
}
