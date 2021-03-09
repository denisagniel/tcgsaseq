#'Plotting raw p-values histogram
#'
#'Display the histogram of raw p-values for diagnostic plots
#'
#'@param x a vector of raw p-values
#'
#'@author Boris Hejblum
#'
#'@return a \code{\link[ggplot2]{ggplot}} object
#'
#'@import ggplot2
#'@export
#'
plot_hist_pvals <- function(pvals, binwidth = 0.02){
  
  stopifnot(is.numeric(pvals))
  
  df2plot <- cbind.data.frame("Rawpvals" = pvals)
  
  ggp <- ggplot(df2plot) + 
    geom_histogram(aes(x = Rawpvals), color="white", binwidth = binwidth) +
    xlab("Raw p-values") +
    theme_bw() +
    ggtitle("Histogram of raw p-values", subtitle = "before multiple testing correction")
  
  return(ggp)
}
