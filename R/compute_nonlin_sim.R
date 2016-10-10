#' Computing nonlinear simulations results
#'
#'@examples
#'\dontrun{
#'for (n in c(50,100,150)) {
#'  betas <- seq(-2, 2, length = 11)
#'  standard.l <- list()
#'  for (i in 1:11) {
#'    standard.sim <- nonlin_sim_fn(type = 'standard',
#'                                  nGenes = 100,
#'                                  n = n,
#'                                  beta = betas[i],
#'                                  n_t = 5,
#'                                  re_sd = 1)
#'    standard.l[[as.character(betas[i])]] <- standard.sim
#'  }
#'}
#'}
#'
#'@keywords internal
#'@importFrom stats sd
#@importFrom edgeR DGEList
#@importFrom limma roast duplicateCorrelation
#'@export


nonlin_sim_fn <- function(n = 250,
                          n_t = 5,
                          beta = 2,
                          nGenes = 1000,
                          rho = 0.5,
                          re_sd = 1,
                          type) {

  if(!requireNamespace("limma", quietly=TRUE)){
    stop("Package 'limma' is not available.\n  -> Try running 'install.packages(\"limma\")'\n")
  }else if(!requireNamespace("edgeR", quietly=TRUE)){
    stop("Package 'edgeR' is not available.\n  -> Try running 'install.packages(\"edgeR\")'\n")
  }else{
    y <- -1

    while(any(y < 0)) {
      sim_data <- sim_nonlin_data(n = n, n_t = n_t,
                                  beta = beta, nGenes = nGenes,
                                  type = type,
                                  rho = rho, re_sd = re_sd)

      set_size <- 10
      set_ind <- 1:set_size
      groups <- floor(nGenes/set_size)

      x <- sim_data$x
      y <- sim_data$y
      tt <- sim_data$tt
      indiv <- sim_data$indiv

      design_r <- cbind(x, tt)
    }

    w <- sp_weights(x = x,
                    y = t(y),
                    phi = matrix(tt, ncol = 1),
                    preprocessed = TRUE,
                    doPlot = FALSE)
    voom_w <- voom_weights(x = design_r,
                           y = t(y),
                           preprocessed = TRUE,
                           doPlot = FALSE)

    cor_limma <- limma::duplicateCorrelation(t(y), design_r, block = indiv)

    genesets <- data.frame(Symbol = rep(1:groups, each = set_size),
                           Chr = rep(1:groups, each = set_size))
    # browser()
    ydge <- edgeR::DGEList(counts=t(y), genes=genesets)
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
    tcgp <- vc_test_perm(y = t(y_set),
                         x = x,
                         indiv = indiv,
                         phi = matrix(tt),
                         Sigma_xi = 1,
                         w = w_set)

    roast_tcg <- limma::roast(y=t(y_set),
                              # set.statistic = 'msq',
                              design=design_r,
                              contrast=3,
                              block = indiv,
                              correlation = cor_limma$consensus.correlation,
                              weights = w_set)
    roast_voom <- limma::roast(y=t(y_set),
                               # set.statistic = 'msq',
                               design=design_r,
                               contrast=3,
                               block = indiv,
                               correlation = cor_limma$consensus.correlation,
                               weights = voomw_set)

    roast_edger <- limma::roast(y=ydge,
                                index = index,
                                # set.statistic = 'msq',
                                design=design_r,
                                contrast=3,
                                correlation = cor_limma$consensus.correlation,
                                block = indiv)
    deseq <- deseq_fn_gs(y, x, tt, indiv, set_ind)
    pvals <- c(
      'tcgsaseq'=tcg$pval,
      'tcgsaseq_perm'=tcgp$pval,
      'roast_tcg'=roast_tcg$p.value["Mixed", "P.Value"],
      'roast_voom'=roast_voom$p.value["Mixed", "P.Value"],
      'roast_edger'=roast_edger$PValue.Mixed,
      'deseq'=deseq
    )
    pvals
  }
}

#'@keywords internal
sim_nonlin_data <- function(n = 250,
                            n_t = 5,
                            beta = 2,
                            nGenes = 1000,
                            rho = 0.5,
                            re_sd = 1,
                            gene_sd = 1,
                            type) {
  x <- cbind(1, matrix(rnorm(n*n_t, mean = 100, sd = 50), n*n_t, 1))
  u <- rexp(nGenes, rate = 1/100)
  e <- t(replicate(nGenes, rnorm(n*n_t)))*u + u
  r.ints <- replicate(n, rnorm(nGenes, sd = 0.01*apply(e, 1, stats::sd)))
  b.0 <- t(apply(r.ints, 1, rep, each = n_t))
  y <- t(e + b.0) + rowMeans(x)
  y <- y*rowMeans(y)/1000
  tt <- matrix(runif(n*n_t), ncol = 1)

  b.1 <- matrix(rep(rnorm(n, sd = re_sd), each = n_t), n*n_t, 1)
  bb.1 <- t(t(b.1) + beta)
  y <- y + rowSums(tt * bb.1)
  y <- ifelse(y <= 0, 1e-7, y)
  indiv <- rep(1:n, each = n_t)

  # lcpm <- apply(y, MARGIN=2,function(v){log2((v+0.5)/(sum(v)+1)*10^6)})
  lcpm <- log2((y+0.5)/(colSums(y) + 1)*10^6)
  lcpm <- pmax(ceiling(lcpm), 1)

  list(x = x, tt = tt, y = lcpm, indiv = indiv)
}
