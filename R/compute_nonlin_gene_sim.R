#' Computing nonlinear simulations results
#'
#'@examples
#'
#'\dontrun{
#'
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
#'
#'}
#'
#'@keywords internal
#@importFrom DESeq2 DESeqDataSetFromMatrix DESeq results estimateDispersionsGeneEst nbinomLRT dispersions
#@importFrom edgeR DGEList calcNormFactors
#@importFrom limma roast duplicateCorrelation
#'@export


nonlin_sim_fn <- function(n = 100,
                          beta = 2,
                          nGenes = 100,
                          rho = 0.5,
                          re_sd = 1,
                          p0 = 0.9) {
  if(!requireNamespace("DESeq2", quietly=TRUE)){
    stop("Package 'DESeq2' is not available.\n  -> Try installing it from Bioconductor\n")
  }else if(!requireNamespace("limma", quietly=TRUE)){
    stop("Package 'limma' is not available.\n  -> Try installing it from Bioconductor\n")
  }else if(!requireNamespace("edgeR", quietly=TRUE)){
    stop("Package 'edgeR' is not available.\n  -> Try installing it from Bioconductor\n")
  }else{

    # browser()
    y <- -1
    while(any(y < 0)) {
      sim_data <- sim_nonlin_data(n = n,
                                  beta = beta,
                                  nGenes = nGenes,
                                  rho = rho,
                                  re_sd = re_sd,
                                  p0 = p0)

      x <- sim_data$x
      y <- sim_data$y
      g <- sim_data$grp
      indiv <- sim_data$indiv

      design_r <- cbind(x, g)
    }
    # browser()
    gene_w <- sp_weights(x = x,
                         y = t(y),
                         phi = matrix(g, ncol = 1),
                         preprocessed = TRUE,
                         doPlot = FALSE,
                         gene_based = TRUE)
    i_w <- sp_weights(x = x,
                      y = t(y),
                      phi = matrix(g, ncol = 1),
                      preprocessed = TRUE,
                      doPlot = FALSE,
                      gene_based = FALSE)
    it_w <- sp_weights(x = x,
                       y = t(y),
                       phi = matrix(g, ncol = 1),
                       preprocessed = TRUE,
                       doPlot = FALSE,
                       gene_based = FALSE,
                       transform = TRUE)
    voom_w <- voom_weights(x = design_r,
                           y = t(y),
                           preprocessed = TRUE,
                           doPlot = FALSE)

    # browser()
    ydge <- edgeR::DGEList(counts=exp(t(y)))
    ydge <- edgeR::calcNormFactors(ydge)
    ydge <- edgeR::estimateDisp(ydge, design_r, robust=TRUE)
    edger_fit <- edgeR::glmQLFit(ydge, design = design_r)
    edger_res <- edgeR::glmQLFTest(edger_fit,coef=3)
    # browser()
    tcgg <- vc_test_asym(y = t(y),
                         x = x,
                         indiv = indiv,
                         phi = matrix(g),
                         Sigma_xi = 1,
                         w = gene_w)
    tcgi <- Dmisc::myTry(vc_test_asym(y = t(y),
                                      x = x,
                                      indiv = indiv,
                                      phi = matrix(g),
                                      Sigma_xi = 1,
                                      w = i_w))
    tcgt <- Dmisc::myTry(vc_test_asym(y = t(y),
                                      x = x,
                                      indiv = indiv,
                                      phi = matrix(g),
                                      Sigma_xi = 1,
                                      w = it_w))


    # tcgg_perm <- vc_test_perm(y = t(y),
    #                           x = x,
    #                           indiv = indiv,
    #                           phi = matrix(g),
    #                           Sigma_xi = 1,
    #                           w = gene_w)
    # tcgi_perm <- vc_test_perm(y = t(y_set),
    #                           x = x,
    #                           indiv = indiv,
    #                           phi = matrix(tt),
    #                           Sigma_xi = 1,
    #                           w = iw_set)
    # tcgt_perm <- vc_test_perm(y = t(y_set),
    #                           x = x,
    #                           indiv = indiv,
    #                           phi = matrix(tt),
    #                           Sigma_xi = 1,
    #                           w = itw_set)

    # browser()
    # index <- list(X=(1:nGenes %in% yinds))
    # tcggp <- tcgip <- tcgtp <- rep(NA, ncol(y))
    # for (i in 1:ncol(y)) {
    #   y_set <- y[,i, drop = FALSE]
    #   gw_set <- gene_w[i,,drop = FALSE]
    #   iw_set <- i_w[i,,drop = FALSE]
    #   itw_set <- it_w[i,,drop = FALSE]
    #   voomw_set <- voom_w[i,,drop = FALSE]
    #   tcgg <- Dmisc::myTry(vc_test_asym(y = t(y_set),
    #                                     x = x,
    #                                     indiv = indiv,
    #                                     phi = matrix(tt),
    #                                     Sigma_xi = 1,
    #                                     w = gw_set))
    #   if (!Dmisc::isErr(tcgg)) tcggp[i] <- tcgg$pval
    #   tcgi <- Dmisc::myTry(vc_test_asym(y = t(y_set),
    #                                     x = x,
    #                                     indiv = indiv,
    #                                     phi = matrix(tt),
    #                                     Sigma_xi = 1,
    #                                     w = iw_set))
    #   if (!Dmisc::isErr(tcgi)) tcgip[i] <- tcgi$pval
    #   tcgt <- Dmisc::myTry(vc_test_asym(y = t(y_set),
    #                                     x = x,
    #                                     indiv = indiv,
    #                                     phi = matrix(tt),
    #                                     Sigma_xi = 1,
    #                                     w = itw_set))
    #   if (!Dmisc::isErr(tcgt)) tcgtp[i] <- tcgt$pval
    #
    # }
    #   # tcggp <- tcgip <- tcgtp <- NA
    #       tcgg <- Dmisc::myTry(vc_test_asym(y = t(y_set),
    #                           x = x,
    #                           indiv = indiv,
    #                           phi = matrix(tt),
    #                           Sigma_xi = 1,
    #                           w = gw_set))
    #       if (!Dmisc::isErr(tcgg)) tcggp[i] <- tcgg$pval
    #     tcgi <- Dmisc::myTry(vc_test_asym(y = t(y_set),
    #                          x = x,
    #                          indiv = indiv,
    #                          phi = matrix(tt),
    #                          Sigma_xi = 1,
    #                          w = iw_set))
    #     if (!Dmisc::isErr(tcgi)) tcgip[i] <- tcgi$pval
    #       tcgt <- Dmisc::myTry(vc_test_asym(y = t(y_set),
    #                           x = x,
    #                           indiv = indiv,
    #                           phi = matrix(tt),
    #                           Sigma_xi = 1,
    #                           w = itw_set))
    #       if (!Dmisc::isErr(tcgt)) tcgtp <- tcgt$pval

    #

    #   roast <- limma::roast
    #   roast_tcg <- limma::roast(y=t(y_set),
    #                              # set.statistic = 'msq',
    #                              design=design_r,
    #                              contrast=3,
    #                              block = indiv,
    #                              correlation = cor_limma$consensus.correlation,
    #                              weights = iw_set)
    #   roast_voom <- limma::roast(y=t(y_set),
    #                              # set.statistic = 'msq',
    #                              design=design_r,
    #                              contrast=3,
    #                              block = indiv,
    #                              correlation = cor_limma$consensus.correlation,
    #                              weights = voomw_set)
    #
    #   roast_edger <- limma::roast(y=ydge,
    #                                       index = index,
    #                                       # set.statistic = 'msq',
    #                                       design=design_r,
    #                                       contrast=3,
    #                               correlation = cor_limma$consensus.correlation,
    #                                       block = indiv)
    # browser()
    vv <- limma::voom(ydge,design_r)
    voom_fit <- limma::lmFit(vv,design_r)
    voom_fit <- limma::eBayes(voom_fit)

    deseq <- deseq_fn(round(exp(y)), x, g, indiv)
    pvals <- data.frame(
      'tcgsaseq_g'=tcgg$gene_pval,
      # 'tcgsaseq_gperm'=tcgg_perm$pval,
      'tcgsaseq_i'=tcgi$gene_pval,
      # 'tcgsaseq_iperm'=tcgi_perm$pval,
      'tcgsaseq_it'=tcgt$gene_pval,
      # 'tcgsaseq_itperm'=tcgt_perm$pval,
      'limma_voom'=voom_fit$p.value[,3],
      # 'roast_tcg'=roast_tcg$p.value["Mixed", "P.Value"],
      # 'roast_voom'=roast_voom$p.value["Mixed", "P.Value"],
      'edger'=edger_res$table[,'PValue'],
      'deseq'=deseq
    )
    pvals
  }

  sim_nonlin_data <- function(n = 250,
                              beta = 2,
                              nGenes = 1000,
                              rho = 0.5,
                              re_sd = 1,
                              p0 = 0.9) {
    # browser()
    x <- cbind(1, matrix(rnorm(n, mean = 100, sd = 50), n, 1))
    u <- rexp(nGenes, rate = 1/100)
    e <- t(replicate(nGenes, rnorm(n)))*u + u
    b.0 <- replicate(n, rnorm(nGenes, sd = 0.01*apply(e, 1, sd)))
    eta <- t(e + b.0) + rowMeans(x)
    eta <- eta*rowMeans(eta)/1000
    grp <- sample(0:1, size = n, replace = TRUE)

    b.1 <- matrix(rnorm(n, sd = re_sd), n, 1)
    bb.1 <- t(t(b.1) + beta)

    n.h1 <- round(nGenes*(1-p0))

    y <- cbind(eta[,1:n.h1] + matrix(grp * bb.1, nrow(eta), n.h1),
               eta[,-(1:n.h1)])
    y <- ifelse(y <= 0, 1e-7, y)
    indiv <- 1:n

    # lcpm <- apply(y, MARGIN=2,function(v){log2((v+0.5)/(sum(v)+1)*10^6)})
    lcpm <- log2((y+0.5)/(colSums(y) + 1)*10^6)
    lcpm <- pmax(ceiling(lcpm), 1)

    list(x = x, grp = grp, y = lcpm, indiv = indiv)
  }


  Kappa_j <- function(r, alpha_DEseq){
    sig <- stats::qnorm(1-alpha_DEseq/2)
    temp_int <- stats::integrate(function(x){exp(-x^2/2)*stats::pnorm((r*x-sig)/sqrt(1-r^2))},
                                 lower=-sig,
                                 upper=sig)
    1/log(1-alpha_DEseq)*log(1-(1/(1-alpha_DEseq)*sqrt(2/pi)*temp_int$value))
  }

  deseq_fn <- function(y, x, tt, indiv, ind) {

    if(!requireNamespace("DESeq2", quietly=TRUE)){
      stop("Package 'DESeq2' is not available.\n  -> Try installing it from Bioconductor\n")
    }else{
      y_dsq <- DESeq2::DESeqDataSetFromMatrix(countData = t(y),
                                              colData = cbind.data.frame("indiv"=as.factor(indiv),
                                                                         "grp"=tt,
                                                                         "x"=as.numeric(x[,2])),
                                              design = ~ x + grp)
      res_dsq <- DESeq2::DESeq(y_dsq, test="LRT", reduced = ~ x)
      pvals <- DESeq2::results(res_dsq)$pvalue
      # library(DESeq2)
      # y_dsq <- DESeq2::estimateSizeFactors(y_dsq)
      # y_dsq <- DESeq2::estimateDispersionsGeneEst(y_dsq)
      # DESeq2::dispersions(y_dsq) <- mcols(y_dsq)$dispGeneEst
      # # res_dsq <- try(DESeq2::DESeq(y_dsq, test="LRT", reduced = ~ x))
      # res_dsq <- DESeq2::nbinomLRT(y_dsq, reduced = ~ x)
      # # if (class(res_dsq) == 'try-error') res_dsq <- try(DESeq2::DESeq(y_dsq, test="LRT", reduced = ~ x, fitType = 'mean'))
      # pvals <- DESeq2::results(res_dsq)$pvalue
      #   if (length(ind) > 1) {
      #     pmin <- min(pvals[ind], na.rm = TRUE)
      #     Rij <- cor(y[,ind])
      #     Rj <- sapply(1:ncol(Rij), function(j){max(abs(Rij[1:(j-1),j]), na.rm = TRUE)})
      #     if(length(which(is.infinite(Rj)))>0){
      #       Rj[which(is.infinite(Rj))] <- 0
      #     }
      #
      #     Keff <- 1 + sum(sapply(Rj[-1], Kappa_j, alpha_DEseq = 0.05))
      #     Keff_approx <- 1 + sum(sqrt(1-Rj[-1]^(-1.31*log10(0.05))))
      #     Meff <- 1 + sum(1-cor(y[,ind])^2, na.rm = TRUE)/ncol(Rij) #Cheverud-Nyholt
      #     cbind("minTest_exact" = 1-(1-pmin)^Keff,
      #           "minTest_approx" = 1-(1-pmin)^Keff_approx,
      #           "minTest_CN" = 1-(1-pmin)^Meff)
      #   }
      return(pvals)
    }
  }
}
