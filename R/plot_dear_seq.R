#'Plot of gene-wise p-values
#'
#'This function prints the sorted exact p-values along with the Benjamini-Hochberg
#'limit and the 5% threshold
#'
#'@param pvals a vector of length \code{n} containing the raw p-values for
#'each gene
#'
#'@return a plot of sorted gene-wise p-values
#'@import viridisLite
#'@export

plot_dear_seq <- function(pvals, ...){
  t <- c(1:length(pvals))
  s <- (t/length(pvals))*0.05
  df_plot_perm <- data.frame("y" = sort(pvals), "x" = c(1:length(pvals)))
  ggplot()+ scale_y_log10()+
    geom_point(data = df_plot_perm,aes(x = x, y = y, color = viridis(4)[1]), size = 0.5)+
    geom_line(data = df_plot_perm, aes(y = s, x = x,color = viridis(4)[2]), size = 0.5) +
    geom_line(data = df_plot_perm, aes(y = 0.05, x = x, color = "red"), size = 0.5) +
    scale_color_manual(name = "", labels = c("B-H limit", "p-values", "5% threshold"),
                       values = c(viridis(4)[c(1,2)], "red")) + xlab("rank") + 
                       ylab("log10 scale") + xlim(0, length(df_plot_perm$y)) +
                       theme_bw()
}
