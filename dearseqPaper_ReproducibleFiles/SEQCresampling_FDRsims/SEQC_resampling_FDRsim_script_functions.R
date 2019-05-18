rm(list=ls())
library(statmod)
library(limma)
library(edgeR)
library(stats4)
library(parallel)
library(BiocGenerics)
library(S4Vectors)
library(IRanges)
library(GenomeInfoDb)
library(GenomicRanges)
library(Biobase)
library(SummarizedExperiment)
library(backports)
library(DESeq2)
library(tcgsaseq)
library(matrixStats)
library(mvtnorm)

load("SEQC_A12345_B12345_counts.RData")
seed <- 1
sampl_sizes2test=300

gen_1sample <- function(orig_counts_mat){

  n <- nrow(orig_counts_mat)

  cov_mat <- cov(t(orig_counts_mat))

  pert <- mvtnorm::rmvnorm(1, mean = rep(0, n),
                  sigma=cov_mat)

  res <- pmax(round(orig_counts_mat[, sample(1, 1:ncol(orig_counts_mat))] + pert), 0)[1,]
  names(res) <- rownames(orig_counts_mat)
  return(res)
}

gen_counts <- function(n_samples=1, orig_counts_mat){
  return(replicate(n_samples, gen_1sample(orig_counts_mat)))
}


sims_fromSEQCA <- function(orig_counts_mat, sample_size = 10){
  stopifnot(sample_size/2 == floor(sample_size/2))
  sims_counts <- gen_counts(n_samples = sample_size, orig_counts_mat)
  group <- factor(c(rep(0, sample_size/2), rep(1, sample_size/2)))

  # edgeR ----
  # biocLite("edgeR")
  y <- DGEList(counts = sims_counts, genes=rownames(sims_counts))
  isexpr <- rowSums(cpm(y)>1) >= 6
  hasannot <- rowSums(is.na(y$genes))==0
  if(sum(isexpr & hasannot)>0){
    y <- y[isexpr & hasannot, , keep.lib.sizes=FALSE]
  }
  y <- calcNormFactors(y)
  Group2 <- group
  design <- model.matrix(~ group)
  y <- estimateDisp(y, design, robust=TRUE)
  fit <- glmQLFit(y, design, robust=TRUE)
  qlf <- glmQLFTest(fit)
  #topTags(qlf,n=15)
  edger_pv <- qlf$table[,4]

  # DESeq2 ----
  names(Group2) <- colnames(y$counts)
  y_dsq <- DESeq2::DESeqDataSetFromMatrix(countData = y$counts,
                                          colData =  cbind.data.frame("group2"=Group2),
                                          design = ~ group2)
  res_dsq <- DESeq2::DESeq(y_dsq, test="LRT", reduced = ~ 1)
  deseq2_pv <- DESeq2::results(res_dsq)$pvalue


  # varseq ----
  ints <- matrix(1,nrow=ncol(sims_counts),ncol=1)
  phi <- as.numeric(as.character(Group2))
  if(sum(isexpr & hasannot)>0){
    sims_counts_dearseq <- sims_counts[isexpr & hasannot, ]
  }else{
    sims_counts_dearseq <- sims_counts
  }

  tcg_pv <- dearseq::varseq(covariates=ints, exprmat=sims_counts_dearseq, variables2test=matrix(phi, ncol=1),
                             gene_based_weights=FALSE, which_weights="loclin",
                             which_test="asymptotic",
                             preprocessed=FALSE, doPlot = FALSE)$pval

  tcgperm_pv <- dearseq::varseq(covariates=ints, exprmat=sims_counts_dearseq, variables2test=matrix(phi, ncol=1),
                                 gene_based_weights=FALSE, which_weights="loclin",
                                 which_test="permutation", n_perm = 1000,
                                 preprocessed=FALSE, doPlot = FALSE)$pval
  # voom:limma ----
  y <- DGEList(counts=sims_counts, genes=rownames(sims_counts))
  if(sum(isexpr & hasannot)>0){
    y <- y[isexpr & hasannot, , keep.lib.sizes=FALSE]
  }
  vv <- limma::voom(y, design, plot=FALSE)
  voom_fit <- limma::lmFit(vv,design)
  voom_fit <- limma::eBayes(voom_fit)


  pvs <- data.frame(edger = edger_pv, deseq2 = deseq2_pv,
                    voom = voom_fit$p.value[, 2],
                    tcg_asym = tcg_pv[, "rawPval"],
                    tcg_perm = tcgperm_pv[, "rawPval"]
  )
  size <- apply(pvs, 2, function(x){sum(x<0.05, na.rm = TRUE)/length(x)})
  fdrs <- apply(apply(pvs, 2, p.adjust, "BH"), 2, function(x){1*(sum(x<0.05, na.rm = TRUE)>0)})

  return(cbind("FDR" = fdrs, "TypeIerr" = size))
}


run_sims <- function(all_orig_counts,
                     ngenes = 2000,
                     sampl_sizes2test = 10
){
  stopifnot(ngenes <= nrow(all_orig_counts))

  res <- list()

  orig_counts <- all_orig_counts[matrixStats::rowSds(all_orig_counts) != 0, ]
  orig_counts_sub <- orig_counts[sample(size=ngenes, x=1:nrow(orig_counts)), ]
  orig_counts_subno0 <- orig_counts_sub[rowSums(orig_counts_sub) != 0, ]
  res_temp <- sims_fromSEQCA(orig_counts_mat = orig_counts_subno0, sample_size = sampl_sizes2test)
  res[["FDR"]] <- res_temp[, "FDR"]
  res[["TypeIerr"]] <- res_temp[, "TypeIerr"]

  return(res)
}

Asamples <- SEQC_A12345_B12345[, grep("A", colnames(SEQC_A12345_B12345))]
set.seed(seed*2017)
test_sims <- run_sims(all_orig_counts=Asamples, sampl_sizes2test = sampl_sizes2test)
test_sims

write.table(t(as.data.frame(test_sims$FDR)), file=paste0("tcgsaseq_FDRsims_sampsize", sampl_sizes2test, "_seed", seed, ".txt"),
            row.names = FALSE, sep="\t")
write.table(t(as.data.frame(test_sims$TypeIerr)), file=paste0("tcgsaseq_TIerrsims_sampsize", sampl_sizes2test, "_seed", seed, ".txt"),
            row.names = FALSE, sep="\t")

