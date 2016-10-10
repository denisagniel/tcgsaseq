#'Time-course Gene Set Analyis
#'
#'Wrapper function for performing gene set analysis of longitudinal RNA-seq data
#'
#'@param y a numeric matrix of size \code{G x n} containing the raw RNA-seq counts or
#'preprocessed expressions from \code{n} samples for \code{G} genes.
#'
#'@param x a numeric matrix of size \code{n x p} containing the model covariates from
#'\code{n} samples (design matrix).
#'
#'@param phi a numeric design matrix of size \code{n x K} containing the \code{K} variables
#'to be tested
#'
#'@param genesets either a vector of index or subscripts that defines which columns of \code{y}
#'constitute the invesigated geneset. Can also be a \code{list} of index when several genesets are
#'tested at once, such as the first element of a \code{\link[GSA:GSA.read.gmt]{gmt}} object.
#'If \code{NULL}, then genewise p-values are returned.
#'
#'@param indiv a vector of length \code{n} containing the information for
#'attributing each sample to one of the studied individuals. Coerced
#'to be a \code{factor}.
#'
#'@param Sigma_xi a matrix of size \code{K x K} containing the covariance matrix
#'of the \code{K} random effects.
#'
#'@param which_weights a character string indicating which method to use to estimate
#'the mean-variance relationship wheights. Possibilities are \code{"loclin"},
#'\code{"voom"} or \code{"none"} (in which case no weighting is performed).
#'Default is \code{"loclin"}.
#'See \code{\link{sp_weights}} and \code{\link{voom_weights}} for details.
#'
#'@param which_test a character string indicating which method to use to approximate
#'the variance component score test, either \code{"permutation"} or \code{"asymptotic"}.
#'Default is \code{"permutation"}.
#'
#'@param n_perm the number of perturbations
#'
#'@param preprocessed a logical flag indicating whether the expression data have
#'already been preprocessed (e.g. log2 transformed). Default is \code{FALSE}, in
#'which case \code{y} is assumed to contain raw counts and is normalized into
#'log(counts) per million.
#'
#'@param doPlot a logical flag indicating whether the mean-variance plot should be drawn.
#' Default is \code{FALSE}.
#'
#'@param gene_based a logical flag used for "loclin" weights, indicating whether to estimate
#'weights at the gene-level. Default is \code{FALSE}, when weights will be estimated at the
#'observation-level.
#'
#'@param bw a character string indicating the smoothing bandwidth selection method to use. See
#'\code{\link[stats]{bandwidth}} for details. Possible values are \code{"ucv"}, \code{"SJ"},
#'\code{"bcv"}, \code{"nrd"} or \code{"nrd0"}
#'
#'@param kernel a character string indicating which kernel should be used.
#'Possibilities are \code{"gaussian"}, \code{"epanechnikov"}, \code{"rectangular"},
#'\code{"triangular"}, \code{"biweight"}, \code{"tricube"}, \code{"cosine"},
#'\code{"optcosine"}. Default is \code{"gaussian"} (NB: \code{"tricube"} kernel
#'corresponds to the loess method).
#'
#'@param exact a logical flag indicating whether the non-parametric weights accounting
#'for the mean-variance relationship should be computed exactly or extrapolated
#'from the interpolation of local regression of the mean against the
#'variance. Default is \code{FALSE}, which uses interporlation (faster computation).
#'
#'@param transform a logical flag used for "loclin" weights, indicating whether values should be
#'transformed to uniform for the purpose of local linear smoothing. This may be helpful if tail
#'observations are sparse and the specified bandwidth gives suboptimal performance there.
#'Default is \code{FALSE}.
#'
#'@param padjust_methods multiple testing correction method used if \code{genesets}
#'is a list. Default is "BH", i.e. Benjamini-Hochberg procedure for contolling the FDR.
#'Other possibilities are: \code{"holm"}, \code{"hochberg"}, \code{"hommel"},
#'\code{"bonferroni"} or \code{"BY"} (for Benjamini-Yekutieli procedure).
#'
#'@param lowess_span smoother span for the lowess function, between 0 and 1. This gives
#'the proportion of points in the plot which influence the smooth at each value.
#'Larger values give more smoothness. Only used if \code{which_weights} is \code{"voom"}.
#'Default is \code{0.5}.
#'
#'@param homogen_traj a logical flag indicating whether trajectories should be considered homogeneous.
#'Default is \code{FALSE} in which case trajectories are not only tested for trend, but also for heterogeneity.
#'
#'@return A list with the following elements:\itemize{
#'   \item \code{which_test}:
#'   \item \code{preprocessed}:
#'   \item \code{n_perm}:
#'   \item \code{pval}: associated p-value
#' }
#'
#'@seealso \code{\link{sp_weights}} \code{\link{vc_test_perm}} \code{\link{vc_test_asym}} \code{\link{p.adjust}}
#'
#'@references Agniel D, Hejblum BP, Variance component score test for
#'time-course gene set analysis of longitudinal RNA-seq data, \emph{submitted}, 2016.
#'
#'@references Law, C. W., Chen, Y., Shi, W., & Smyth, G. K. (2014). voom: Precision
#'weights unlock linear model analysis tools for RNA-seq read counts. \emph{Genome
#'Biology}, 15(2), R29.
#'
#'@importFrom stats p.adjust
#'
#'@examples
#'#rm(list=ls())
#'n <- 200
#'r <- 12
#'t <- matrix(rep(1:3), 4, ncol=1, nrow=r)
#'sigma <- 0.4
#'b0 <- 1
#'
#'#under the null:
#'b1 <- 0
#'
#'y.tilde <- b0 + b1*t + rnorm(r, sd = sigma)
#'y <- t(matrix(rnorm(n*r, sd = sqrt(sigma*abs(y.tilde))), ncol=n, nrow=r) +
#'       matrix(rep(y.tilde, n), ncol=n, nrow=r))
#'x <- matrix(1, ncol=1, nrow=r)
#'
#'#run test
#'res <- tcgsa_seq(y, x, phi=t, genesets=lapply(0:9, function(x){x*10+(1:10)}),
#'                         Sigma_xi=matrix(1), indiv=rep(1:4, each=3), which_test="asymptotic",
#'                         which_weights="none", preprocessed=TRUE)
#'res_genes <- tcgsa_seq(y, x, phi=t, genesets=NULL,
#'                       Sigma_xi=matrix(1), indiv=rep(1:4, each=3), which_test="asymptotic",
#'                       which_weights="none", preprocessed=TRUE)
#'quantile(res_genes$pvals[, "rawPval"])
#'@export
tcgsa_seq <- function(y, x, phi, genesets,
                      indiv = rep(1, nrow(x)), Sigma_xi = diag(ncol(phi)),
                      which_test = c("permutation", "asymptotic"),
                      which_weights = c("loclin", "voom", "none"),
                      n_perm = 1000,
                      preprocessed = FALSE, doPlot = TRUE, gene_based = FALSE,
                      bw = "nrd",
                      kernel = c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "tricube", "cosine", "optcosine"),
                      exact = FALSE, transform = FALSE,
                      padjust_methods = c("BH", "BY", "holm", "hochberg", "hommel", "bonferroni"),
                      lowess_span = 0.5,
                      homogen_traj = FALSE){


  stopifnot(is.matrix(y))

  if(!preprocessed){
    y_lcpm <- t(apply(y, MARGIN=1, function(v){log2((v+0.5)/(sum(v)+1)*10^6)}))
    preprocessed <- TRUE
  }else{
    y_lcpm <- y
  }
  rm(y)

  if(is.data.frame(x)){
    warning("design matrix 'x' is a data.frame instead of a matrix:\n all variables (including factors) are converted to numeric...")
    x <- as.matrix(as.data.frame(lapply(x, as.numeric)))
  }

  if(is.data.frame(phi)){
    warning("design matrix 'phi' is a data.frame instead of a matrix:\n all variables (including factors) are converted to numeric... ")
    phi <- as.matrix(as.data.frame(lapply(phi, as.numeric)))
  }

  if(det(crossprod(cbind(x, phi)))==0){
    stop("crossprod(x, phi) cannot be inversed. x and phi are likely colinear...")
  }

  if(length(padjust_methods)>1){
    padjust_methods <- padjust_methods[1]
  }
  stopifnot(padjust_methods %in% c("BH", "BY", "holm", "hochberg", "hommel", "bonferroni"))

  if(length(which_weights)>1){
    which_weights <- which_weights[1]
  }
  stopifnot(which_weights %in% c("loclin", "voom", "none"))

  if(length(which_test)>1){
    which_test <- which_test[1]
  }
  stopifnot(which_test %in% c("asymptotic", "permutation"))

  w <-  switch(which_weights,
               loclin = sp_weights(y = y_lcpm, x = x, phi=phi,
                                   preprocessed = preprocessed, doPlot=doPlot,
                                   gene_based = gene_based,
                                   bw = bw, kernel = kernel,
                                   exact = exact, transform = transform),
               voom = voom_weights(y = y_lcpm, x = cbind(x, phi),
                                   preprocessed = preprocessed, doPlot = doPlot,
                                   lowess_span = lowess_span),
               none = matrix(1, ncol=ncol(y_lcpm), nrow=nrow(y_lcpm))
  )


  if(which_test == "asymptotic"){
    n_perm <- NA
  }

  if(is.null(genesets)){
    cat("'genesets' argument not provided => only genewise p-values are computed\n")

    if(which_test == "asymptotic"){
      rawPvals <- vc_test_asym(y = y_lcpm, x = x, indiv = indiv, phi = phi,
                               w = w, Sigma_xi = Sigma_xi,
                               genewise_pvals = TRUE, homogen_traj = homogen_traj)$gene_pvals
    }else if(which_test == "permutation"){
      rawPvals <- vc_test_perm(y = y_lcpm, x = x, indiv = indiv, phi = phi,
                               w = w, Sigma_xi = Sigma_xi,
                               n_perm=n_perm, genewise_pvals = TRUE, homogen_traj = homogen_traj)$gene_pvals
    }

    pvals <- data.frame("rawPval" = rawPvals, "adjPval" = stats::p.adjust(rawPvals, padjust_methods))
    if(!is.null(rownames(y_lcpm))){
      rownames(pvals) <- rownames(y_lcpm)
    }
  }else if(is.list(genesets)){
    if(class(genesets[[1]])=="character"){
      gene_names_measured <- rownames(y_lcpm)
      prop_meas <- sapply(genesets, function(x){length(intersect(x, gene_names_measured))/length(x)})
      if(sum(prop_meas)!=length(prop_meas)){
        warning("Some genes in the investigated genesets were not measured:\nremoving those genes from the geneset definition...")
        genesets <- lapply(genesets, function(x){x[which(x %in% gene_names_measured)]})
      }
    }else if(!is.vector(genesets)){
      stop("'genesets' argument provided but is neither a list nor an atomic vector")
    }

    if(which_test == "asymptotic"){

        rawPvals <- sapply(genesets, FUN = function(gs){
          vc_test_asym(y = y_lcpm[gs, ], x = x, indiv = indiv, phi = phi,
                       w = w[gs, ], Sigma_xi = Sigma_xi,
                       genewise_pvals = FALSE, homogen_traj = homogen_traj)$set_pval}
        )
    } else if(which_test == "permutation"){
        rawPvals <- sapply(genesets, FUN = function(gs){
          vc_test_perm(y = y_lcpm[gs, ], x = x, indiv = indiv, phi = phi,
                       w = w[gs, ], Sigma_xi = Sigma_xi,
                       n_perm=n_perm, genewise_pvals = FALSE, homogen_traj = homogen_traj)$set_pval}
        )
      }

    pvals <- data.frame("rawPval" = rawPvals, "adjPval" = stats::p.adjust(rawPvals, padjust_methods))
    if(!is.null(names(genesets))){
      rownames(pvals) <- names(genesets)
    }

  }else{

    if(!is.vector(genesets)){
      stop("'genesets' argument provided but is neither a list nor a vector")
    }

    if(class(genesets)=="character"){
      gene_names_measured <- rownames(y_lcpm)
      if((length(intersect(genesets, gene_names_measured))/length(x)) != 1){
        warning("Some genes in the investigated genesets were not measured:\n removing those genes from the geneset definition...")
        genesets <- genesets[which(genesets %in% gene_names_measured)]
      }
    }

    res_test <- switch(which_test,
                       asymptotic = vc_test_asym(y = y_lcpm[genesets, ], x = x, indiv = indiv, phi = phi,
                                                 w = w[genesets, ], Sigma_xi = Sigma_xi,
                                                 genewise_pvals = FALSE, homogen_traj = homogen_traj),
                       permutation = vc_test_perm(y = y_lcpm[genesets, ], x = x, indiv = indiv, phi = phi,
                                                  w = w[genesets, ], Sigma_xi = Sigma_xi, n_perm = n_perm,
                                                  genewise_pvals = FALSE, homogen_traj = homogen_traj)
    )
    pvals <- data.frame("rawPval" = res_test$set_pval, "adjPval" = NA)
    padjust_methods <- NA

  }

  return(list("which_test" = which_test, "preprocessed" = preprocessed, "n_perm" = n_perm,
              "genesets" = genesets, "pvals" = pvals))

}
