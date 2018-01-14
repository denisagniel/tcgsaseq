#' Computing negative binomial simulations results
#'
#'@examples
#'\dontrun{
#'for (n in c(50,100, 150)){
#'standard.l <- list()
#'betas <- seq(-2, 2, length = 11)
#'stdl.l <- list()
#'ind <- 0
#'for (i in 1:11){
#'  ind <- ind + 1
#'
#'  if (betas[i] == 0){
#'    gsd <- rsd <- 0
#'  }else{
#'    gsd <- rsd <- 1
#'  }
#'
#'  standard.sim <- nb_sim_fn(type = 'standard',
#'                            nGenes = 100,
#'                            n = n,
#'                            beta = betas[i],
#'                            re_sd = rsd,
#'                            gene_sd = gsd,
#'                            n_t = 5)
#'
#'  stdl.l[[ind]] <- c(standard.sim, beta = betas[i])
#'}
#'standard.l <- do.call(rbind, stdl.l)
#'}
#'}
#'
#'
#'@keywords internal
#'@importFrom stats rnorm rexp rnbinom qnorm integrate pnorm
#@importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#@importFrom edgeR DGEList estimateDisp calcNormFactors results DESeq
#@importFrom limma roast
#'@export

nb_sim_fn <- function(n = 250,
                      n_t = 5,
                      beta = 2,
                      nGenes = 1000,
                      re_sd = 1,
                      gene_sd = 1,
                      type){

  if(!requireNamespace("DESeq2", quietly=TRUE)){
    stop("Package 'DESeq2' is not available.\n  -> Try installing it from Bioconductor\n")
  }else if(!requireNamespace("limma", quietly=TRUE)){
    stop("Package 'limma' is not available.\n  -> Try installing it from Bioconductor\n")
  }else if(!requireNamespace("edgeR", quietly=TRUE)){
    stop("Package 'edgeR' is not available.\n  -> Try installing it from Bioconductor\n")
  }else{

    sim_data <- sim_nb_data(n = n, n_t = n_t, re_sd = re_sd,
                            beta = beta, nGenes = nGenes, type = type,
                            gene_sd = gene_sd)

    set_size <- 10
    set_ind <- 1:set_size
    groups <- floor(nGenes/set_size)

    x <- sim_data$x
    y <- sim_data$y
    tt <- sim_data$tt
    indiv <- sim_data$indiv

    design_r <- cbind(x, tt)

    w <- sp_weights(x = x,
                    y = t(y),
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
    w_set <- abs(w[set_ind,])
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
                              correlation = 0,
                              weights = w_set)
    roast_voom <- limma::roast(y=t(y_set),
                               # set.statistic = 'msq',
                               design=design_r,
                               contrast=3,
                               block = indiv,
                               correlation = 0,
                               weights = voomw_set)
    roast_edger <- limma::roast(y=ydge,
                                index = index,
                                # set.statistic = 'msq',
                                design=design_r,
                                contrast=3,
                                block = indiv,
                                correlation = 0)
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
sim_nb_data <- function(n = 250,
                        n_t = 5,
                        beta = 2,
                        nGenes = 1000,
                        re_sd = 1,
                        gene_sd = 1,
                        type) {
  stopifnot(type %in% c('standard' , 'td_variance'))
  mu_x <- rexp(n)*10
  # mu_x <- rnorm(n)
  x <- cbind(1, rep(rnorm(n, mean = mu_x), each = n_t))
  tt <- c(replicate(n, sample(1:10, size = n_t)))
  # tt <- c(replicate(n, 1:n_t))
  b_0 <- replicate(nGenes, rep(rnorm(n), each = n_t))
  b_1 <- replicate(nGenes, rep(rnorm(n, sd = re_sd), each = n_t))
  gene_f <- rnorm(nGenes)*gene_sd
  mu <- y <- b_0
  int <- 1000
  for (ii in 1:ncol(b_0)) {
    mu[,ii] <- int + b_0[,ii] + rowSums(x) + (b_1[,ii] + beta + gene_f[ii])*(tt*x[,2])
    mu[,ii] <- ifelse(mu[,ii] < 0, 0, mu[,ii])
    if (type == 'standard') {
      y[,ii] <- rnbinom(n*n_t, mu = mu[,ii], size = 1) + 1
    } else if (type == 'td_variance') {
      y[,ii] <- rnbinom(n*n_t, mu = mu[,ii], size = 1) + 1
    }
  }
  indiv <- rep(1:n, each = n_t)
  list(x = x, tt = tt, y = y, y.0 = y[mu != 0], indiv = indiv)
}


#'@keywords internal
#'
Kappa_j <- function(r, alpha_DEseq){
  sig <- stats::qnorm(1-alpha_DEseq/2)
  temp_int <- stats::integrate(function(x){exp(-x^2/2)*stats::pnorm((r*x-sig)/sqrt(1-r^2))},
                               lower=-sig,
                               upper=sig)
  1/log(1-alpha_DEseq)*log(1-(1/(1-alpha_DEseq)*sqrt(2/pi)*temp_int$value))
}


#'@keywords internal
deseq_fn_gs <- function(y, x, tt, indiv, ind) {

  if(!requireNamespace("DESeq2", quietly=TRUE)){
    stop("Package 'DESeq2' is not available.\n  -> Try installing it from Bioconductor\n")
  }else{
    y_dsq <- DESeq2::DESeqDataSetFromMatrix(countData = t(y),
                                            colData = cbind.data.frame("indiv"=as.factor(indiv),
                                                                       "time"=as.numeric(tt),
                                                                       "x"=as.numeric(x[,2])),
                                            design = ~ x + time)
    res_dsq <- DESeq2::DESeq(y_dsq, test="LRT", reduced = ~ x)
    pvals <- DESeq2::results(res_dsq)$pvalue

    pmin <- min(pvals[ind], na.rm = TRUE)
    Rij <- cor(y[,ind])
    Rj <- sapply(1:ncol(Rij), function(j){max(abs(Rij[1:(j-1),j]), na.rm = TRUE)})
    if(length(which(is.infinite(Rj)))>0){
      Rj[which(is.infinite(Rj))] <- 0
    }

    Keff <- 1 + sum(sapply(Rj[-1], Kappa_j, alpha_DEseq = 0.05))
    Keff_approx <- 1 + sum(sqrt(1-Rj[-1]^(-1.31*log10(0.05))))
    Meff <- 1 + sum(1-cor(y[,ind])^2, na.rm = TRUE)/ncol(Rij) #Cheverud-Nyholt
    cbind("minTest_exact" = 1-(1-pmin)^Keff,
          "minTest_approx" = 1-(1-pmin)^Keff_approx,
          "minTest_CN" = 1-(1-pmin)^Meff)
  }
}
