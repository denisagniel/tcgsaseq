#' Computing simulations results observation-wise
#'
#'@examples
#'
#'\dontrun{
#'rm(list=ls())
#'data_sims <- data_sim_voomlike(seed=1, do_gs=FALSE)
#'res <- compute_sim_voomlike_genewise(counts = data_sims$counts,
#'                                     design = data_sims$design,
#'                                     indiv = data_sims$indiv,
#'                                     alternative=TRUE,
#'                                     fixed_eff = 0.5,
#'                                     fixed_eff_sd = 0,
#'                                     rand_eff_sd = 0,
#'                                     RE_indiv_sd=NULL)
#'res_all <- cbind(res$res_voom, res$res_perso, res$res_noweights, res$res_DEseq, res$res_edgeR)
#'colnames(res_all) <- c(paste0(rep(c("asym", "perm", "camera", "roast"), 3),
#'                              rep(c("_voom", "_perso", "_noweights"), each=4)),
#'                       paste0("DESeq2_minTest_", c("exact", "approx", "CN")), "roast_edgeR")
#'}
#'
#'@keywords internal
#'@importFrom stats rnorm qnorm pnorm integrate
#@importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#@importFrom edgeR DGEList estimateDisp calcNormFactors roast.DGEList
#@importFrom limma roast camera
#'@export
compute_sim_voomlike_genewise <- function(counts, design, indiv, alternative=FALSE,
                                          fixed_eff = 0.5, fixed_eff_sd = 0, #0.2,
                                          rand_eff_sd = 0.25, RE_indiv_sd=NULL, eps_sd=0.05){

  if(!requireNamespace("DESeq2", quietly=TRUE)){
    stop("Package 'DESeq2' is not available.\n  -> Try installing it from Bioconductor\n")
  }else if(!requireNamespace("limma", quietly=TRUE)){
    stop("Package 'limma' is not available.\n  -> Try installing it from Bioconductor\n")
  }else if(!requireNamespace("edgeR", quietly=TRUE)){
    stop("Package 'edgeR' is not available.\n  -> Try installing it from Bioconductor\n")
  }else{

    logcoutspm <- apply(counts, MARGIN=2, function(v){log2((v + 0.5)/(sum(v) + 1)*10^6)})

    if(alternative){
      n <- ncol(counts)
      nb_g <- nrow(counts)
      beta_time <- stats::rnorm(1, mean=fixed_eff, sd=fixed_eff_sd)
      gamma_genes_time <- matrix(stats::rnorm(nb_g, mean=0, sd=rand_eff_sd), ncol=1)
      logcoutspm_alt <- logcoutspm + (gamma_genes_time + beta_time)%*%design[, "time"]

      # for(gs in gs_keep){
      #   nb_g_gs <- length(gs)
      #   gamma_genes_time <- matrix(stats::rnorm(nb_g_gs, mean = 0, sd = rand_eff_sd), ncol = 1)
      #   logcoutspm_alt[gs,] <- logcoutspm[gs,] + (gamma_genes_time+beta_time)%*%design[, "time"]
      # }

      err <- matrix(stats::rnorm(n*nb_g, mean = 0, sd = eps_sd), nrow = nb_g, ncol = n)
      logcoutspm_alt <- logcoutspm_alt + err

    }else{
      logcoutspm_alt <- logcoutspm
    }

    #Individual Random effects
    if(!is.null(RE_indiv_sd)){
      RE_indiv <- stats::rnorm(n=length(unique(indiv)), mean=0, sd=RE_indiv_sd)
      logcoutspm_alt <- logcoutspm_alt + RE_indiv[as.numeric(as.character(indiv))]
    }

    # Weights estimation ----
    #########################

    gene_w <- sp_weights(x = design[,-which(colnames(design)=="time")],
                         y = logcoutspm_alt,
                         phi =  cbind(design[, "time"]),
                         preprocessed = TRUE,
                         doPlot = FALSE,
                         gene_based = TRUE)
    i_w <- sp_weights(x = design[,-which(colnames(design)=="time")],
                      y = logcoutspm_alt,
                      phi =  cbind(design[, "time"]),
                      preprocessed = TRUE,
                      doPlot = FALSE,
                      gene_based = FALSE)
    it_w <- sp_weights(x = design[,-which(colnames(design)=="time")],
                       y = logcoutspm_alt,
                       phi =  cbind(design[, "time"]),
                       preprocessed = TRUE,
                       doPlot = FALSE,
                       gene_based = FALSE,
                       transform = TRUE)
    voom_w <- voom_weights(x = design,
                           y = logcoutspm_alt,
                           preprocessed = TRUE,
                           doPlot = FALSE)


    # edgeR ----
    ###########
    ydge <- edgeR::DGEList(counts=floor(exp(logcoutspm_alt)), genes=rownames(logcoutspm_alt))
    ydge <- edgeR::calcNormFactors(ydge)
    ydge <- edgeR::estimateDisp(ydge, design, robust=TRUE)
    edger_fit <- edgeR::glmQLFit(ydge, design = design)
    edger_res <- edgeR::glmQLFTest(edger_fit, coef=which(colnames(design)=="time"))


    # # Perm ----
    # ###########
    # yinds <- 1:10
    # y_set <- logcoutspm_alt[yinds,, drop = FALSE]
    # gw_set <- gene_w[yinds, ,drop = FALSE]
    # iw_set <- i_w[yinds, ,drop = FALSE]
    # itw_set <- it_w[yinds, ,drop = FALSE]
    # voomw_set <- voom_w[yinds, ,drop = FALSE]
    #
    #
    # tcgg_perm <- vc_test_perm(y = y_set,
    #                           x = design[,-which(colnames(design)=="time")],
    #                           indiv = indiv,
    #                           phi = cbind(design[, "time"]),
    #                           Sigma_xi = 1,
    #                           w = gw_set)
    # tcgi_perm <- vc_test_perm(y = y_set,
    #                           x = design[,-which(colnames(design)=="time")],
    #                           indiv = indiv,
    #                           phi = cbind(design[, "time"]),
    #                           Sigma_xi = 1,
    #                           w = iw_set)
    # tcgt_perm <- vc_test_perm(y = y_set,
    #                           x = design[,-which(colnames(design)=="time")],
    #                           indiv = indiv,
    #                           phi = cbind(design[, "time"]),
    #                           Sigma_xi = 1,
    #                           w = itw_set)


    # tcgsaseq ----
    ############
    tcggp <- tcgip <- tcgtp <- rep(NA, nb_g)
    tcgg_permp <- tcgi_permp <- tcgt_permp <- rep(NA, nb_g)
    for (i in 1:nb_g) {
      y_set <- logcoutspm_alt[i,, drop = FALSE]
      gw_set <- gene_w[i,,drop = FALSE]
      iw_set <- i_w[i,,drop = FALSE]
      itw_set <- it_w[i,,drop = FALSE]
      voomw_set <- voom_w[i,,drop = FALSE]
      tcgg <- try(vc_test_asym(y = y_set,
                               x = design[,-which(colnames(design)=="time")],
                               indiv = indiv,
                               phi = cbind(design[, "time"]),
                               Sigma_xi = 1,
                               w = gw_set))
      if(!inherits(tcgg, "try-error")){tcggp[i] <- tcgg$pval}
      tcgg_perm <- try(vc_test_perm(y = y_set,
                                x = design[,-which(colnames(design)=="time")],
                                indiv = indiv,
                                phi = cbind(design[, "time"]),
                                Sigma_xi = 1,
                                w = gw_set))
      if(!inherits(tcgg_perm, "try-error")){tcgg_permp[i] <- tcgg_perm$pval}

      tcgi <- try(vc_test_asym(y = y_set,
                               x = design[,-which(colnames(design)=="time")],
                               indiv = indiv,
                               phi = cbind(design[, "time"]),
                               Sigma_xi = 1,
                               w = iw_set))
      if(!inherits(tcgi, "try-error")){tcgip[i] <- tcgi$pval}
      tcgi_perm <- try(vc_test_perm(y = y_set,
                                    x = design[,-which(colnames(design)=="time")],
                                    indiv = indiv,
                                    phi = cbind(design[, "time"]),
                                    Sigma_xi = 1,
                                    w = iw_set))
      if(!inherits(tcgi_perm, "try-error")){tcgi_permp[i] <- tcgi_perm$pval}

      tcgt <- try(vc_test_asym(y = y_set,
                               x = design[,-which(colnames(design)=="time")],
                               indiv = indiv,
                               phi = cbind(design[, "time"]),
                               Sigma_xi = 1,
                               w = itw_set))
      if(!inherits(tcgt, "try-error")){tcgtp[i] <- tcgt$pval}
      tcgt_perm <- try(vc_test_perm(y = y_set,
                                    x = design[,-which(colnames(design)=="time")],
                                    indiv = indiv,
                                    phi = cbind(design[, "time"]),
                                    Sigma_xi = 1,
                                    w = itw_set))
      if(!inherits(tcgt_perm, "try-error")){tcgt_permp[i] <- tcgt_perm$pval}
    }

    # voom ----
    ############
    vv <- limma::voom(ydge, design)
    voom_fit <- limma::lmFit(vv, design)
    voom_fit <- limma::eBayes(voom_fit)


    # DESeq2 ----
    #############
    deseq <- deseq_fn(y = floor(exp(logcoutspm_alt)), x = design[,-which(colnames(design) %in% c("(Intercept)", "time")), drop=FALSE],
                      phi = cbind(design[, "time"]), indiv = indiv)


    # wrapping results
    pvals <- data.frame(
      'tcgsaseq_g'=tcggp,
      'tcgsaseq_gperm'=tcgg_perm$pval,
      'tcgsaseq_i'=tcgip,
      'tcgsaseq_iperm'=tcgi_perm$pval,
      'tcgsaseq_it'=tcgtp,
      'tcgsaseq_itperm'=tcgt_perm$pval,
      'limma_voom'=voom_fit$p.value[,3],
      'edger'=edger_res$table[,'PValue'],
      'deseq'=deseq
    )
    return(pvals)

  }
}
