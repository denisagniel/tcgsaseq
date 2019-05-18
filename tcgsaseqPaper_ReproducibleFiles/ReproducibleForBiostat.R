# Kidney Transplant study ####
# PBT gene sets
# Dorr et al. 2015
###################
rm(list=ls())

# Data import
# data available at : http://conservancy.umn.edu/handle/11299/177631.1
FPKMs <- as.matrix(read.delim("http://conservancy.umn.edu/bitstream/handle/11299/177631.1/FPKMs.txt", row.names=1))
clinical <- read.delim("http://conservancy.umn.edu/bitstream/handle/11299/177631.1/clinical.txt")

# loading gene sets of interest
data(PBTreproducible)

# data preprocess
FPKMs <- FPKMs[, as.character(clinical$ID)]
cs_fpkm <- colSums(FPKMs)
TPMs <- (FPKMs/cs_fpkm)*10^6
design <- as.matrix(clinical[, -1])
log2TPMs <- log2(TPMs + 1)
data("PBTrepoducible_gmt")

# limma preprocessing
library(limma)
gs_ind_list <- limma::ids2indices(PBT_gmt[[1]], identifier=rownames(log2TPMs), remove.empty=FALSE)
cor_limma <- limma::duplicateCorrelation(log2TPMs, design[, -which(colnames(design)=="SID")],
                                         block = design[,"SID"])

# weights computation
w <- tcgsaseq::sp_weights(x=design[,-which(colnames(design)=="month")],
                          y=log2TPMs,
                          phi = cbind(design[, "month"]),
                          doPlot = TRUE,
                          exact=FALSE,
                          preprocessed = TRUE,
                          gene_based=TRUE
)
w_voom <- tcgsaseq::voom_weights(y=log2TPMs, x=design, doPlot = TRUE, preprocessed = TRUE)


set.seed(54321)
library(doParallel)
cl <- parallel::makeCluster(parallel::detectCores()-2, outfile="")
doParallel::registerDoParallel(cl)
res_all_par <- foreach(gs=1:length(PBT_gmt[[1]]), .packages=c("tcgsaseq", "limma")) %dopar% {

  cat(gs, "/", length(PBT_gmt[[1]]), sep="")
  gs_test <- intersect(PBT_gmt[[1]][[gs]], rownames(log2TPMs))

  if(length(gs_test)<2){
    res_voomRoast <- NA
    res_tcgsaseq <- NA
    warning("Less than 2 genes measured: output is NA")
  }else{
    y_test <- log2TPMs[gs_test,]
    x_test <- design[,-which(colnames(design)%in%c("month", "SID"))]
    w_test_perso <- w[gs_test, ]
    w_test_voom <- w_voom[gs_test, ]
    phi_test <- cbind(design[, "month"])
    indiv_test <- cbind(design[, "SID"])

    n <- ncol(y_test)
    g <- length(gs_test)
    y_T_vect <- as.vector(y_test)
    indiv_vect <- as.factor(rep(design[, "SID"], g))
    g_vect <- as.factor(rep(gs_test, n))
    x_vect <-do.call(rbind, replicate(g, x_test, simplify=FALSE))
    phi_vect <- do.call(rbind, replicate(g, phi_test, simplify=FALSE))

    res_voomRoast <- limma::roast(y=log2TPMs, index=gs_ind_list[[gs]],
                                  design=design[, -2], contrast=1,
                                  block = factor(indiv_test), correlation = cor_limma$consensus.correlation,
                                  weights = w_voom)$p.value["Mixed", "P.Value"]

    res_tcgsaseq <- vc_test_perm(y = y_test, x=x_test, indiv=indiv_test, phi=phi_test,
                                 w = w_test_perso, n_perm=5000,
                                 Sigma_xi = as.matrix(diag(ncol(phi_test))))[["set_pval"]]

  }
  cat(" DONE\n")
  return(list("voomRoast" = res_voomRoast, "tcgsaseq" = res_tcgsaseq))
}

stopCluster(cl)

res_all_par_df <- NULL
for(nam in names(res_all_par[[1]])){
  res_all_par_df <- cbind(res_all_par_df, sapply(res_all_par, "[[", nam))
}
colnames(res_all_par_df) <- names(res_all_par[[1]])


# Results  & plots
###################
#rm(list=ls())
#load("PBT_gmt.RData")
#load("res_all_parList_kidney_PBT_54321.RData")
library(reshape2)
library(ggplot2)

sapply(lapply(res_all_par, function(x){x<0.05}), colMeans, na.rm=TRUE)
PBT_gmt$geneset.description[which(res_all_par_df[,"tcgsaseq"]<0.05)]
PBT_gmt$geneset.description[which(res_all_par_df[, "voomRoast"]<0.05)]


res_adj <- apply(res_all_par_df, 2, FUN=p.adjust, method="BH")
colMeans(apply(res_adj, 2, function(x){x<0.1}), na.rm=TRUE)
PBT_gmt$geneset.description[which(res_adj[, "tcgsaseq"]<0.1)]
PBT_gmt$geneset.description[which(res_adj[, "voomRoast"]<0.1)]

res2plot <- cbind.data.frame("Geneset" = unlist(PBT_gmt$geneset.description),
                             "tcgsaseq"=res_all_par_df[, "tcgsaseq"],
                             "voom+ROAST"=res_all_par_df[, "voomRoast"])
res2plot$Geneset <- factor(res2plot$Geneset, levels=res2plot$Geneset, ordered = TRUE)
levels(res2plot$Geneset) <- gsub("Quantitative ", "", gsub("NK", "Natural Killer", gsub("CTL", "Cytotoxic T cells", gsub("B-cell", "B cell", gsub("cell", "cells", gsub(" Induced", "", gsub(fixed=TRUE, " (DSA) selective", "", gsub("New ", "", gsub(" associated", "", gsub("-Associated", "", gsub(" transcript", "", gsub(" transcripts", "", gsub(" Transcripts", "", gsub(" burden", "", levels(res2plot$Geneset)))))))))))))))
res2plot_melted <- reshape2::melt(res2plot, id.vars = "Geneset", variable.name = "Method", value.name = "Pval")
p <- (ggplot(res2plot_melted)
      + geom_bar(aes(x=Geneset, y=Pval, group=Method, color=Method, fill=Method), stat="identity", position = "dodge", alpha=0.7)
      + coord_flip()
      + ylim(0,0.6)
      + scale_fill_manual("\nMethod", values=c("#d7191c", "#00004d"), breaks=c("voom+ROAST", "tcgsaseq"))
      + scale_color_manual("\nMethod", values=c("#d7191c", "#00004d"), breaks=c("voom+ROAST", "tcgsaseq"))
      + theme_bw()
      + xlab("PBT Gene set")
      + ylab("p-value")
      + theme(legend.key = element_rect(colour = "white"),
              panel.grid.major.y=element_blank(),
              axis.title=element_text(size=23),
              axis.text=element_text(size=17),
              legend.title=element_text(size=20),
              legend.text=element_text(size = 16),
              legend.key.size=unit(0.32,"inches"))
      + scale_linetype_manual(name="Type I error", values = c(2))
      + geom_hline(aes(yintercept=0.05, linetype="5%"), size=0.65,
                   col="black", show.legend = TRUE)
      + guides(linetype=guide_legend("\nSignificance\nthreshold"),
               fill=guide_legend("\nMethod", override.aes = list(linetype=0)))
)
#pdf(file="KidneyPBT_pval_barchart.pdf", width=12, height=6)
p
#dev.off()






# Baduel et al. Study -----

rm(list=ls())
library(tcgsaseq)

##import the data
#################

data(baduel_5gs)

set.seed(54321)
KAvsTBG <- tcgsa_seq(y=log2(expr_norm_corr+1),
                     x=apply(as.matrix(design[, c("Intercept", "Vernalized", "AgeWeeks",
                                                  "Vernalized_Population", "AgeWeeks_Population"),
                                              drop=FALSE]),2, as.numeric),
                     phi=as.matrix(design[, c("PopulationKA"), drop=FALSE]),
                     genesets=baduel_gmt$genesets,
                     which_test = "permutation", which_weights = "loclin",
                     n_perm=200, preprocessed = TRUE, doPlot = TRUE, gene_based=TRUE)
KAvsTBG$pvals$rawPval

set.seed(54321)
Cold <- tcgsa_seq(y=log2(expr_norm_corr+1),
                  x=apply(as.matrix(design[, c("Intercept", "PopulationKA", "AgeWeeks",
                                               "AgeWeeks_Population"), drop=FALSE]), 2, as.numeric),
                  phi=as.matrix(design[, c("Vernalized", "Vernalized_Population")]),
                  genesets=baduel_gmt$genesets,
                  which_test = "permutation", which_weights = "loclin",
                  n_perm=200, preprocessed = TRUE, doPlot = TRUE, gene_based=TRUE)
Cold$pvals$rawPval
