#' Computing simulations results
#'
#'@examples
#'
#'\dontrun{
#'data_sims <- data_sim_voomlike(seed=1)
#'res <- compute_sim_voomlike(counts = data_sims$counts,
#'                            design = data_sims$design,
#'                            gs_keep = data_sims$gs_keep,
#'                            indiv = data_sims$indiv,
#'                            alternative=TRUE,
#'                            fixed_eff = 0.5,
#'                            fixed_eff_sd = 0,
#'                            rand_eff_sd = 0,
#'                            RE_indiv_sd=NULL)
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
compute_sim_voomlike <- function(counts, design, gs_keep, indiv, alternative=FALSE,
                                 fixed_eff = 0.5, fixed_eff_sd = 0, #0.2,
                                 rand_eff_sd = 0.25, RE_indiv_sd=NULL, eps_sd=0.05, alpha_DEseq=0.05){

  if(!requireNamespace("DESeq2", quietly=TRUE)){
    stop("Package 'DESeq2' is not available.\n  -> Try running 'install.packages(\"DESeq2\")'\n")
  }else if(!requireNamespace("limma", quietly=TRUE)){
    stop("Package 'limma' is not available.\n  -> Try running 'install.packages(\"limma\")'\n")
  }else if(!requireNamespace("edgeR", quietly=TRUE)){
    stop("Package 'edgeR' is not available.\n  -> Try running 'install.packages(\"edgeR\")'\n")
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

    #png(width=5.5, height=4.5, units="in", file="voom_ex.png", res=300)
    #w_voom <- voom_weights(x=design, y=t(counts2), doPlot = TRUE, preprocessed = FALSE)
    #dev.off()
    w_voom <- voom_weights(x=design, y=logcoutspm_alt, doPlot = FALSE, preprocessed = TRUE)

    w_perso <- sp_weights(x=design[,-which(colnames(design)=="time")],
                          y=logcoutspm_alt,
                          phi = cbind(design[, "time"]),
                          doPlot = FALSE,
                          exact=FALSE,
                          preprocessed = TRUE
    )

    # DESeq ----
    ############
    DESeq_data <- DESeq2::DESeqDataSetFromMatrix(countData = floor(exp(logcoutspm_alt)),
                                                 colData = cbind.data.frame("indiv"=as.factor(indiv),
                                                                            "time"=as.numeric(design[, "time"]),
                                                                            "group"=as.factor(design[, "group2"])),
                                                 design = ~ group + time)
    DEseq_res <- DESeq2::DESeq(DESeq_data, test="LRT", reduced = ~ group)
    DEseq_univ_res_pvals <- DESeq2::results(DEseq_res)$pvalue
    names(DEseq_univ_res_pvals) <- rownames(logcoutspm_alt)

    Kappa_j <- function(r, alpha_DEseq){
      sig <- stats::qnorm(1-alpha_DEseq/2)
      temp_int <- stats::integrate(function(x){exp(-x^2/2)*stats::pnorm((r*x-sig)/sqrt(1-r^2))},
                            lower=-sig,
                            upper=sig)
      1/log(1-alpha_DEseq)*log(1-(1/(1-alpha_DEseq)*sqrt(2/pi)*temp_int$value))
    }



    # edgeR ----
    ############
    edgeR_data <- edgeR::DGEList(counts=floor(exp(logcoutspm_alt)), genes = rownames(logcoutspm_alt))
    edgeR_data <- edgeR::calcNormFactors(edgeR_data)
    edgeR_data <- edgeR::estimateDisp(edgeR_data,
                                      design=cbind("time"=as.numeric(design[, "time"]),
                                                   "group"=as.numeric(design[, "group2"])),
                                      robust=TRUE)
    #roast <- limma::roast

    # Testing ----
    ##############
    res_voom <- NULL
    res_perso <- NULL
    res_noweights <- NULL
    res_DEseq <- NULL
    res_edgeR <- NULL
    gs_ind_list <- limma::ids2indices(gs_keep, identifier=rownames(logcoutspm_alt))
    cor_limma <- limma::duplicateCorrelation(logcoutspm_alt, design, block = indiv)

    for(gs in 1:length(gs_keep)){
      gs_test <- gs_keep[[gs]]
      y_test <- logcoutspm_alt[gs_test,]
      x_test <- design[,-which(colnames(design)=="time")]
      w_test_perso <- w_perso[as.numeric(gs_test), ]
      w_test_voom <- w_voom[as.numeric(gs_test), ]
      phi_test <- cbind(design[, "time"])

      n <- ncol(y_test)
      g <- length(gs_test)
      y_T_vect <- as.vector(y_test)
      indiv_vect <- rep(indiv, g)
      g_vect <- as.factor(rep(1:g, n))
      x_vect <-do.call(rbind, replicate(g, x_test, simplify=FALSE))
      phi_vect <- do.call(rbind, replicate(g, phi_test, simplify=FALSE))

      y_deseq <- floor(exp(y_test))
      DEseq_pmin <- min(DEseq_univ_res_pvals[gs_test], na.rm = TRUE)
      Rij <- cor(t(y_deseq))
      Rj <- sapply(1:g, function(j){max(abs(Rij[1:(j-1),j]), na.rm = TRUE)})
      if(length(which(is.infinite(Rj)))>0){
        Rj[which(is.infinite(Rj))] <- 0
      }

      Keff <- 1 + sum(sapply(Rj[-1], Kappa_j, alpha_DEseq = 0.05))
      Keff_approx <- 1 + sum(sqrt(1-Rj[-1]^(-1.31*log10(alpha_DEseq))))
      Meff <- 1 + sum(1-cor(t(y_deseq))^2, na.rm = TRUE)/g #Cheverud-Nyholt
      res_DEseq <- rbind(res_DEseq, cbind("minTest_exact" = 1-(1-DEseq_pmin)^Keff,
                                          "minTest_approx" = 1-(1-DEseq_pmin)^Keff_approx,
                                          "minTest_CN" = 1-(1-DEseq_pmin)^Meff)
      )

      res_edgeR <- rbind(res_edgeR,
                         edgeR::roast.DGEList(y=edgeR_data, index = gs_ind_list[[gs]],
                                              design = design, contrast = 3)$p.value["Mixed", "P.Value"]
      )
      #fry(y=edgeR_data, index = gs_ind_list[[gs]], design = design, contrast = 3, sort="mixed")

      res_voom <- rbind(res_voom, cbind("asym" = vc_test_asym(y = y_test, x=x_test, indiv=indiv, phi=phi_test,
                                                              Sigma_xi = as.matrix(diag(ncol(phi_test))),
                                                              w = w_test_voom)[["pval"]],
                                        "permut" = vc_test_perm(y = y_test, x=x_test, indiv=indiv, phi=phi_test,
                                                                Sigma_xi = as.matrix(diag(ncol(phi_test))),
                                                                w = w_test_voom)[["pval"]],
                                        "camera" = limma::camera(y=logcoutspm_alt, index=gs_ind_list[[gs]],
                                                                 design=design, contrast=3,
                                                                 weights = w_voom)$PValue,
                                        "roast" = limma::roast(y=logcoutspm_alt, index=gs_ind_list[[gs]],
                                                               design=design, contrast=3,
                                                               block = indiv, correlation = cor_limma$consensus.correlation,
                                                               weights = w_voom)$p.value["Mixed", "P.Value"])
      )
      res_perso <- rbind(res_perso, cbind("asym" = vc_test_asym(y = y_test, x=x_test, indiv=indiv, phi=phi_test,
                                                                Sigma_xi = as.matrix(diag(ncol(phi_test))),
                                                                w = w_test_perso)[["pval"]],
                                          "permut" = vc_test_perm(y = y_test, x=x_test, indiv=indiv, phi=phi_test,
                                                                  Sigma_xi = as.matrix(diag(ncol(phi_test))),
                                                                  w = w_test_perso)[["pval"]],
                                          "camera" = limma::camera(y=logcoutspm_alt, index=gs_ind_list[[gs]],
                                                                   design=design, contrast=3,
                                                                   weights = w_perso)$PValue,
                                          "roast" = limma::roast(y=logcoutspm_alt, index=gs_ind_list[[gs]],
                                                                 design=design, contrast=3,
                                                                 block = indiv, correlation = cor_limma$consensus.correlation,
                                                                 weights = w_perso)$p.value["Mixed", "P.Value"])
      )
      res_noweights <- rbind(res_noweights, cbind("asym" = vc_test_asym(y = y_test, x=x_test, indiv=indiv, phi=phi_test,
                                                                        Sigma_xi = as.matrix(diag(ncol(phi_test))),
                                                                        w = rep(1, length(gs_test)))[["pval"]],
                                                  "permut" = vc_test_perm(y = y_test, x=x_test, indiv=indiv, phi=phi_test,
                                                                          Sigma_xi = as.matrix(diag(ncol(phi_test))),
                                                                          w = rep(1, length(gs_test)))[["pval"]],
                                                  "camera" = limma::camera(y=logcoutspm_alt, index=gs_ind_list[[gs]],
                                                                           design=design, contrast=3)$PValue,
                                                  "roast" = limma::roast(y=logcoutspm_alt, index=gs_ind_list[[gs]],
                                                                         block = indiv, correlation = cor_limma$consensus.correlation,
                                                                         design=design, contrast=3)$p.value["Mixed", "P.Value"])
      )

    }

    return(list("res_voom"=res_voom, "res_perso"=res_perso, "res_noweights"=res_noweights,
                "res_DEseq"=res_DEseq, "res_edgeR"=res_edgeR))
  }
}
