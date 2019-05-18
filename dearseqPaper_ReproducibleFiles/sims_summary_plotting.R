rm(list=ls())
library(ggplot2)

#### Data MGMT ----
res <- readRDS("NB_NonLin_sims_12072018.rds")
load("aggregatedResults_permpe_12122018.RData")

res2plot <- reshape2::melt(res, id.vars=c("setting", "n", "method"), variable.name = "indicator")
res2plot$indicator <- factor(as.character(res2plot$indicator), levels = c("ti_err", "fdr", "pwr", "tpr"), ordered=TRUE)
levels(res2plot$indicator) <- c("Type-I error", "FDR", "Power", "TPR")
res2plot$setting <- factor(as.character(res2plot$setting), levels = c("nb", "nonlin", "seqc_rs"), ordered = TRUE)

temp <- res_agg_permpe
levels(temp$method) <- c("voom", "edger", "dsq", "vs", "vsp")
res2plot <- rbind.data.frame(res2plot, temp[, colnames(res2plot)])

setting_fullname <- c(
  nonlin = "Nonlinear",
  nb = "Negative Binomial",
  seqc_rs = "SEQC resampling + Gaussian noise"
)
indicator_fullname <- c(
  "Type-I error" = "Type-I error",
  Power = "Statistical power",
  TPR = "True Discovery Rate",
  FDR = "False Discovery Rate"
)

#### Output Plots -----

ggplot(dplyr::filter(res2plot, indicator %in% c("Type-I error", "FDR")),
       aes(x=n, y=value)) +
  theme_bw() +
  xlab("Sample size") +
  ylab("Indicator value") +
  ggtitle("Monte-Carlo estimation over 1,000 simulation runs",
          subtitle = "(nominal testing level at 5%)") +
  geom_hline(aes(yintercept=0.05, size="5%"), color="black", linetype=3) +
  geom_point(aes(col=method, shape=method), size=2) +
  geom_line(aes(col=method)) +
  scale_x_continuous(breaks=c(4, 8, 16, 50, 100, 150, 200, 300)) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(size=4.5)) +
  theme(axis.text.y = element_text(size=6)) +
  scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1), minor_breaks = seq(0.125, 0.875, by=0.25), limits = c(0,1)) +
  scale_size_manual("Nominal level", values=c(0.5)) +
  scale_color_manual("Method", labels = c("limma-voom","edgeR", "DESeq2", "dearseq - asymp", "dearseq - perm"),
                     values = c(viridis::viridis(n=4)[-4], "gold2", heat.colors(n=5)[3])) +
  scale_shape("Method", labels = c("limma-voom","edgeR", "DESeq2", "dearseq - asymp", "dearseq - perm")) +

  facet_grid(setting ~ indicator, labeller = labeller(setting = setting_fullname, indicator = indicator_fullname)) +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.box = "horizontal") +
  NULL
ggsave(filename = "FDRTypeI.pdf", width = 7, height = 10)
ggsave(filename = "FDRTypeI.png", width = 7, height = 10, device = "png")



ggplot(dplyr::filter(res2plot, indicator %in% c("Power", "TPR")),
       aes(x=n, y=value)) +
  theme_bw() +
  xlab("Sample size") +
  ylab("Indicator value") +
  ggtitle("Monte-Carlo estimation over 1,000 simulation runs",
          subtitle = "(nominal testing level at 5%)") +
  geom_hline(aes(yintercept=0.05, size="5%"), color="black", linetype=3) +
  geom_point(aes(col=method, shape=method), size=2) +
  geom_line(aes(col=method)) +
  scale_x_continuous(breaks=c(4, 8, 16, 50, 100, 150, 200, 300)) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(size=4.5)) +
  theme(axis.text.y = element_text(size=6)) +
  scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1), minor_breaks = seq(0.125, 0.875, by=0.25), limits = c(0,1)) +
  scale_size_manual("Nominal level", values=c(0.5)) +
  scale_color_manual("Method", labels = c("limma-voom","edgeR", "DESeq2", "dearseq - asymp", "dearseq - perm"),
                     values = c(viridis::viridis(n=4)[-4], "gold2", grDevices::heat.colors(n=5)[3])) +
  scale_shape("Method", labels = c("limma-voom","edgeR", "DESeq2", "dearseq - asymp", "dearseq - perm")) +

  facet_grid(setting ~ indicator, labeller = labeller(setting = setting_fullname, indicator = indicator_fullname)) +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.box = "horizontal") +
  NULL

ggsave(filename = "TPRPower.pdf", width = 7, height = 8)
ggsave(filename = "TPRPower.png", width = 7, height = 8, device = "png")




