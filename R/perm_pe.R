#'Exact permutation p-values
#'
#'Calculates exact p-values for permutation tests when permutations are randomly drawn with replacement. This implementation is based on
#'
#'@author Belinda Phipson and Gordon Smyth (adapted by Boris Hejblum)
#'
#'@param nperm_supobs number of permutations that yielded test statistics at least as extreme as the observed data. Can be a vector or an array of values.
#'@param nperm_eff number of permutations effectively computed.
#'@param total_possible_nperm total number of permutations possible.
#'
#'@references Phipson B, and Smyth GK (2010). Permutation p-values should never be zero:
#'calculating exact p-values when permutations are randomly drawn. \emph{Statistical Applications
#'in Genetics and Molecular Biology}, Volume 9, Issue 1, Article 39.
#'\url{http://www.statsci.org/smyth/pubs/PermPValuesPreprint.pdf}
#'
#'@importFrom statmod gauss.quad.prob
#'@importFrom stats pbinom
#'
#'@export

perm_pe <- function(nperm_supobs, nperm_eff, total_possible_nperm){

  # (slightly) adapted from the implementation by Belinda Phipson and Gordon Smyth from the R package statmod

  nbnodes_quad <- 128
  z <- statmod::gauss.quad.prob(n = nbnodes_quad, l = 0, u = 0.5/total_possible_nperm)
  npvals <- length(nperm_supobs)
  prob <- rep(z$nodes, npvals)
  x2 <- rep(nperm_supobs, each = nbnodes_quad)
  Y <- matrix(data = stats::pbinom(x2, prob = prob, size = nperm_eff), nrow = nbnodes_quad, ncol = npvals)
  int <- 0.5/total_possible_nperm * colSums(z$weights * Y)

  pe <- (nperm_supobs + 1)/(nperm_eff + 1) - int

  return(pe)
}
