#' Small portion of RNA-seq data from plant physiology study.
#'
#'A subsample of the RNA-seq data from Baduel \emph{et al.} studying Arabidopsis Arenosa physiology.
#'
#'@name baduel_5gs
#'@rdname baduel_5gs
#'@aliases baduel baduel_5gs baduel_gmt design expr_norm_corr
#'
#'@usage data(baduel_5gs)
#'
#'@references Baduel P, Arnold B, Weisman CM, Hunter B & Bomblies K, Habitat-Associated Life
#'History and Stress-Tolerance Variation in Arabidopsis Arenosa. \emph{Plant Physiology}: 01875, 2015.
#'
#'@references Agniel D, Hejblum BP, Variance component score test for
#' time-course gene set analysis of longitudinal RNA-seq data, \emph{submitted}, 2016.
#'
#'@format 3 objects\itemize{
#'\item{\code{design}:} a design matrix for the 48 measured samples, containing the following variables:\itemize{
#'  \item \code{SampleName} corresponding column names from \code{expr_norm_corr}
#'  \item \code{Intercept} an intercept variable
#'  \item \code{Population} a factor identifying the plant population
#'  \item \code{Age_weeks} numeric age of the plant at sampling time (in weeks)
#'  \item \code{Replicate} a purely technical variable as replicates are not from the same individual over weeks.
#'  Should not be used in analysis.
#'  \item \code{Vernalized} a logical variable indicating whether the plant had undergone
#'  vernalization (exposition to cold and short day photoperiods)
#'  \item \code{Vernalized} a binary variable indicating whether the plant belonged to the KA
#'  population
#'  \item \code{AgeWeeks_Population} interaction variable between the \code{AgeWeeks} and
#'  \code{Population} variables
#'  \item \code{AgeWeeks_Vernalized} interaction variable between the \code{AgeWeeks} and
#'  \code{Vernalized} variables
#'  \item \code{Vernalized_Population} interaction variable between the \code{Vernalized} and
#'  \code{Population} variables
#'  \item \code{AgeWeeks_Vernalized_Population} interaction variable between the \code{AgeWeeks},
#'  \code{Vernalized} and \code{Population} variables
#'}
#'\item{\code{baduel_gmt}:} a \code{gmt} object containing 5 gene sets of interest (see \code{\link[GSA:GSA.read.gmt]{GSA.read.gmt}})
#'\item{\code{expr_norm_corr}:} a numeric matrix containing the normalized batch corrected expression
#'for the 2454 genes included in either of the 5 gene sets of interests
#'}
#'
#'@examples
#' \dontrun{
#' rm(list=ls())
#' data("baduel_5gs")
#'
#' set.seed(54321)
#' KAvsTBG <- tcgsa_seq(y=log2(expr_norm_corr+1), x=apply(as.matrix(design[, c("Intercept",
#'    "Vernalized", "Age_weeks", "Vernalized_Population", "AgeWeeks_Population"), drop=FALSE]),
#'        2, as.numeric),
#'                      phi=as.matrix(design[, c("PopulationKA"), drop=FALSE]),
#'                      genesets=baduel_gmt$genesets[c(3,5)],
#'                      which_test = "permutation", which_weights = "loclin",
#'                      n_perm=1000, preprocessed = TRUE, doPlot = TRUE)
#'
#' set.seed(54321)
#' Cold <- tcgsa_seq(y=log2(expr_norm_corr+1), x=apply(as.matrix(design[, c("Intercept",
#'    "Age_weeks", "PopulationKA", "AgeWeeks_Population"), drop=FALSE]), 2, as.numeric),
#'                  phi=as.matrix(design[, c("Vernalized", "Vernalized_Population")]),
#'                  genesets=baduel_gmt$genesets[c(3,5)],
#'                  which_test = "permutation", which_weights = "loclin",
#'                  n_perm=1000, preprocessed = TRUE, doPlot = TRUE)
#' }
#'
#'
#' @source \url{http://www.ncbi.nlm.nih.gov/bioproject/PRJNA312410}
#'
#' @keywords datasets
#' @docType data
NULL

