#'Time-course Gene Set Analysis
#'
#'Wrapper function for performing gene set analysis of (potentially longitudinal) RNA-seq data
#'
#'@param y a numeric matrix of size \code{G x n} containing the raw RNA-seq counts or
#'preprocessed expressions from \code{n} samples for \code{G} genes.
#'
#'@param x a numeric matrix of size \code{n x p} containing the model covariates from
#'\code{n} samples (design matrix). Usually, its first column is the intercept (full of
#'\code{1}s).
#'
#'@param phi a numeric design matrix of size \code{n x K} containing the \code{K} variables
#'to be tested
#'
#'@param weights_phi_condi a logical flag indicating whether heteroscedasticity
#'weights computation should be conditional on both the variable(s) to be tested
#'\code{phi} and on covariate(s) \code{x}, or on \code{x} alone. #'Default is \code{TRUE}
#'in which case conditional means are estimated conditionally on both \code{x} and \code{phi}.
#'
#'@param genesets either a vector of index or subscripts that defines which columns of \code{y}
#'constitute the investigated gene set (when only 1 gene set is being tested).
#'Can also be a \code{list} of index (or \code{rownames} of \code{y}) when several
#'gene sets are tested at once, such as the first element of a
#'\code{\link[GSA:GSA.read.gmt]{gmt}} object. If \code{NULL}, then gene-wise p-values are returned.
#'
#'@param indiv a vector of length \code{n} containing the information for
#'attributing each sample to one of the studied individuals. Coerced
#'to be a \code{factor}. Default is \code{NULL} in which case each sample is considered
#'as coming from independent subjects.
#'
#'@param Sigma_xi a matrix of size \code{K x K} containing the covariance matrix
#'of the \code{K} random effects. Only used if \code{homogen_traj} is \code{FALSE}.
#'Default assume diagonal correlation matrix, i.e. independence of random effects.
#'
#'@param which_weights a character string indicating which method to use to estimate
#'the mean-variance relationship weights. Possibilities are \code{"loclin"},
#'\code{"voom"} or \code{"none"} (in which case no weighting is performed).
#'Default is \code{"loclin"}.
#'See \code{\link{sp_weights}} and \code{\link{voom_weights}} for details.
#'
#'@param which_test a character string indicating which method to use to approximate
#'the variance component score test, either \code{"permutation"} or \code{"asymptotic"}.
#'Default is \code{"permutation"}.
#'
#'@param n_perm the number of perturbations. Default is \code{1000}.
#'
#'@param preprocessed a logical flag indicating whether the expression data have
#'already been preprocessed (e.g. log2 transformed). Default is \code{FALSE}, in
#'which case \code{y} is assumed to contain raw counts and is normalized into
#'log(counts) per million.
#'
#'@param doPlot a logical flag indicating whether the mean-variance plot should be drawn.
#' Default is \code{FALSE}.
#'
#'@param gene_based_weights a logical flag used for \code{"loclin"} weights, indicating whether to estimate
#'weights at the gene-level, or rather at the observation-level. Default is \code{TRUE},
#'and weights are then estimated at the gene-level.
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
#'variance. Default is \code{FALSE}, which uses interpolation (faster computation).
#'
#'@param transform a logical flag used for \code{"loclin"} weights, indicating whether values should be
#'transformed to uniform for the purpose of local linear smoothing. This may be helpful if tail
#'observations are sparse and the specified bandwidth gives suboptimal performance there.
#'Default is \code{FALSE}.
#'
#'@param padjust_methods multiple testing correction method used if \code{genesets}
#'is a list. Default is "BH", i.e. Benjamini-Hochberg procedure for controlling the FDR.
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
#'@param na.rm_tcgsaseq logical: should missing values in \code{y} (including
#'\code{NA} and \code{NaN}) be omitted from the calculations?
#'Default is \code{TRUE}.
#'
#'@param verbose logical: should informative messages be printed during the
#'computation? Default is \code{TRUE}.
#'
#'@return A list with the following elements:\itemize{
#'   \item \code{which_test}: a character string carrying forward the value of the '\code{which_test}' argument
#'    indicating which test was perform (either "asymptotic" or "permutation").
#'   \item \code{preprocessed}: a logical flag carrying forward the value of the '\code{preprocessed}' argument
#'   indicating whether the expression data were already preprocessed, or were provided as raw counts and
#'   transformed into log-counts per million.
#'   \item \code{n_perm}: an integer carrying forward the value of the '\code{n_perm}' argument indicating
#'   the number of perturbations performed (\code{NA} if asymptotic test was performed).
#'   \item \code{genesets}: carrying forward the value of the '\code{genesets}' argument defining the gene sets
#'   of interest (\code{NULL} for gene-wise testing).
#'   \item \code{pval}: computed p-values. A \code{data.frame} with one raw for each each gene set, or
#'   for each gene if \code{genesets} argument is \code{NULL}, and with 2 columns: the first one '\code{rawPval}'
#'   contains the raw p-values, the second one contains the FDR adjusted p-values and is either named
#'   '\code{adjPval}' (according to the '\code{padjust_methods}' argument) in the \code{asymptotic} case
#'   or '\code{FDR}' in the \code{permutation} case.
#' }
#'
#'@seealso \code{\link{sp_weights}} \code{\link{vc_test_perm}} \code{\link{vc_test_asym}} \code{\link{p.adjust}}
#'
#'@references Agniel D & Hejblum BP (2017). Variance component score test for
#'time-course gene set analysis of longitudinal RNA-seq data, \emph{Biostatistics},
#'18(4):589-604. \href{https://doi.org/10.1093/biostatistics/kxx005}{10.1093/biostatistics/kxx005}.
#'\href{https://arxiv.org/abs/1605.02351}{arXiv:1605.02351}.
#'
#'@references Law, C. W., Chen, Y., Shi, W., & Smyth, G. K. (2014). voom: Precision
#'weights unlock linear model analysis tools for RNA-seq read counts. \emph{Genome
#'Biology}, 15(2), R29.
#'
#'@importFrom stats p.adjust
#'
#'@examples
#'#rm(list=ls())
#'nsims <- 2 #100
#'res_quant <- list()
#'for(i in 1:2){
#'n <- 2000#0
#'nr <- 3
#'r <- nr*20#4*nr#100*nr
#'t <- matrix(rep(1:nr), r/nr, ncol=1, nrow=r)
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
#'                         Sigma_xi=matrix(1), indiv=rep(1:(r/nr), each=nr), which_test="asymptotic",
#'                         which_weights="none", preprocessed=TRUE)
#'res_genes <- tcgsa_seq(y, x, phi=cbind(t),#, rnorm(r)), #t^2
#'                       genesets=NULL,
#'                       Sigma_xi=diag(1), indiv=rep(1:(r/nr), each=nr), which_test="asymptotic",
#'                       which_weights="none", preprocessed=TRUE)
#'length(res_genes$pvals[, "rawPval"])
#'quantile(res_genes$pvals[, "rawPval"])
#'res_quant[[i]] <- res_genes$pvals[, "rawPval"]
#'}
#'#round(rowMeans(sapply(res_quant, quantile)), 3)
#'#plot(density(unlist(res_quant)))
#'#mean(unlist(res_quant)<0.05)
#'
#'\dontrun{
#'res_genes <- tcgsa_seq(y, x, phi=t, genesets=NULL,
#'                       Sigma_xi=matrix(1), indiv=rep(1:(r/nr), each=nr), which_test="permutation",
#'                       which_weights="none", preprocessed=TRUE, n_perm=1000)
#'
#'mean(res_genes$pvals$rawPval < 0.05)
#'summary(res_genes$pvals$FDR)
#'}
#'@export
tcgsa_seq <- function(y, x, phi, weights_phi_condi = TRUE,
                      genesets,
                      indiv = NULL,
                      Sigma_xi = diag(ncol(phi)),
                      which_test = c("permutation", "asymptotic"),
                      which_weights = c("loclin", "voom", "none"),
                      n_perm = 1000,
                      preprocessed = FALSE, doPlot = TRUE, gene_based_weights = TRUE,
                      bw = "nrd",
                      kernel = c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "tricube", "cosine", "optcosine"),
                      exact = FALSE, transform = FALSE,
                      padjust_methods = c("BH", "BY", "holm", "hochberg", "hommel", "bonferroni"),
                      lowess_span = 0.5,
                      homogen_traj = FALSE,
                      na.rm_tcgsaseq = TRUE,
                      verbose = TRUE){

  stopifnot(is.matrix(y))
  stopifnot(is.matrix(x) | is.data.frame(x))
  stopifnot(is.matrix(phi) | is.data.frame(phi))

  if(sum(is.na(y))>1 & na.rm_tcgsaseq){
    warning(paste("\n\n!!!!!\n'y' contains", sum(is.na(y)), "NA values.",
                  "\nCurrently they are ignored in the computations but you should think carefully about where do those NA/NaN come from...",
                  "\nIf you don't want to ignore those NA/NaN values, set the 'na.rm_tcgsaseq' argument to 'FALSE' (this may lead to errors).\n!!!!!\n"))
  }

  if(!preprocessed){
    y_lcpm <- apply(y, MARGIN=2, function(v){log2((v+0.5)/(sum(v)+1)*10^6)})
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
    stop("crossprod(cbind(x, phi)) cannot be inversed. 'x' and 'phi' are likely colinear...")
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


  if(which_test == "asymptotic"){
    n_perm <- NA

    if(nrow(x) < 10)
      warning("Less than 10 samples: asymptotics likely not reached \nYou should probably run permutation test instead...")
  }



  if(which_test == "permutation"){
    if(is.null(indiv)){
      options(warn = -1)
      N_possible_perms <- factorial(ncol(y_lcpm))
      options(warn = 0)
    }else{
      options(warn = -1)
      N_possible_perms <- prod(sapply(table(indiv), factorial))
      options(warn = 0)
    }

    if(n_perm > N_possible_perms){
      stop(paste("The number of permutations requested 'n_perm' is larger than the total number of existing permutations", N_possible_perms,
                 ". Try a lower number for 'n_perm'"))
    }
  }



  # Computing the weights
  if(which_weights != "none" & verbose){message("Computing the weights... ")}
  w <-  switch(which_weights,
               loclin = sp_weights(y = y_lcpm, x = x, phi = phi, use_phi = weights_phi_condi,
                                   preprocessed = TRUE, doPlot = doPlot,
                                   gene_based = gene_based_weights,
                                   bw = bw, kernel = kernel,
                                   exact = exact, transform = transform, verbose = verbose, na.rm = na.rm_tcgsaseq),
               voom = voom_weights(y = y_lcpm, x = if(weights_phi_condi){cbind(x, phi)}else{x},
                                   preprocessed = TRUE, doPlot = doPlot,
                                   lowess_span = lowess_span),
               none = matrix(1, ncol=ncol(y_lcpm), nrow=nrow(y_lcpm),
                             dimnames = list(rownames(y_lcpm),
                                             colnames(y_lcpm)
                             )
               )
  )
  if(which_weights != "none" & verbose){message("Done!\n")}



  if(is.null(genesets)){
    if(verbose){
      message("'genesets' argument not provided => only gene-wise p-values are computed\n")
    }
    if(which_test == "asymptotic"){
      if(is.null(indiv)){
        indiv <- 1:nrow(x)
      }

      rawPvals <- vc_test_asym(y = y_lcpm, x = x, indiv = indiv, phi = phi,
                               w = w, Sigma_xi = Sigma_xi,
                               genewise_pvals = TRUE, homogen_traj = homogen_traj,
                               na.rm = na.rm_tcgsaseq)$gene_pvals
    }else if(which_test == "permutation"){
      if(is.null(indiv)){
        indiv <- rep(1, nrow(x))
      }

      # constructing residuals
      if(na.rm_tcgsaseq){
        y_lcpm0 <- y_lcpm
        y_lcpm0[is.na(y_lcpm0)] <- 0
        y_lcpm_res <- y_lcpm - t(x%*%solve(crossprod(x))%*%t(x)%*%t(y_lcpm0))
      }else{
        y_lcpm_res <- y_lcpm - t(x%*%solve(crossprod(x))%*%t(x)%*%t(y_lcpm))
      }
      rm(y_lcpm0)
      x_res <- matrix(1, nrow=nrow(x), ncol=1)

      res <- vc_test_perm(y = y_lcpm_res, x = x_res, indiv = indiv, phi = phi,
                          w = w, Sigma_xi = Sigma_xi,
                          n_perm=n_perm, genewise_pvals = TRUE, homogen_traj = homogen_traj,
                          na.rm = na.rm_tcgsaseq)
      rawPvals <- res$gene_pvals

      ds_fdr <- res$fdr
    }

    if (which_test == "permutation"){
      pvals <- data.frame("rawPval" = rawPvals, "FDR"= ds_fdr)
    }
    else if (which_test == "asymptotic"){
      pvals <- data.frame("rawPval" = rawPvals, "adjPval" = stats::p.adjust(rawPvals, padjust_methods))
    }

    if(!is.null(rownames(y_lcpm))){
      rownames(pvals) <- rownames(y_lcpm)
    }
  }else if(is.list(genesets)){

    if(class(genesets[[1]])=="character"){
      if(is.null(rownames(y_lcpm))){
        stop("Gene sets specified as character but no rownames available for the expression matrix")
      }
      gene_names_measured <- rownames(y_lcpm)
      prop_meas <- sapply(genesets, function(x){length(intersect(x, gene_names_measured))/length(x)})
      if(sum(prop_meas)!=length(prop_meas)){
        warning("Some transcripts in the investigated gene sets were not measured:\nremoving those transcripts from the gene set definition...")
        genesets <- lapply(genesets, function(x){x[which(x %in% gene_names_measured)]})
      }
    }

    if(which_test == "asymptotic"){
      if(is.null(indiv)){
        indiv <- 1:nrow(x)
      }

      rawPvals <- sapply(seq_along(genesets), FUN = function(i_gs){
        gs <- genesets[[i_gs]]
        e <- try(y_lcpm[gs, 1, drop=FALSE], silent = TRUE)
        if(inherits(e, "try-error") | length(e)==0){
          warning(paste("Gene set", i_gs, "contains 0 measured transcript: associated p-value cannot be computed"))
          NA
        }else{
          vc_test_asym(y = y_lcpm[gs, , drop=FALSE], x = x, indiv = indiv, phi = phi,
                       w = w[gs, , drop=FALSE], Sigma_xi = Sigma_xi,
                       genewise_pvals = FALSE, homogen_traj = homogen_traj,
                       na.rm = na.rm_tcgsaseq)$set_pval
        }
      }
      )
    } else if(which_test == "permutation"){
      if(is.null(indiv)){
        indiv <- rep(1, nrow(x))
      }

      if(na.rm_tcgsaseq){
        y_lcpm0 <- y_lcpm
        y_lcpm0[is.na(y_lcpm0)] <- 0
        y_lcpm_res <- y_lcpm - t(x%*%solve(crossprod(x))%*%t(x)%*%t(y_lcpm0))
      }else{
        y_lcpm_res <- y_lcpm - t(x%*%solve(crossprod(x))%*%t(x)%*%t(y_lcpm))
      }
      rm(y_lcpm0)

      x_res <- matrix(1, nrow=nrow(x), ncol=1)
      rawPvals <- sapply(seq_along(genesets), FUN = function(i_gs){
        gs <- genesets[[i_gs]]
        e <- try(y_lcpm[gs, 1], silent = TRUE)
        if(inherits(e, "try-error")){
          warning(paste("Gene set", i_gs, "contains 0 measured transcript: associated p-value cannot be computed"))
          NA
        }else{
          vc_test_perm(y = y_lcpm[gs, ], x = x, indiv = indiv, phi = phi,
                       w = w[gs, , drop=FALSE], Sigma_xi = Sigma_xi,
                       n_perm=n_perm, genewise_pvals = FALSE, homogen_traj = homogen_traj,
                       na.rm = na.rm_tcgsaseq)$set_pval
        }
      }
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
      if(is.null(rownames(y_lcpm))){
        stop("Gene sets specified as character but no rownames available for the expression matrix")
      }
      gene_names_measured <- rownames(y_lcpm)
      if((length(intersect(genesets, gene_names_measured))/length(x)) != 1){
        warning("Some transcripts in the investigated gene sets were not measured:\n removing those transcripts from the gene set definition...")
        genesets <- genesets[which(genesets %in% gene_names_measured)]
      }
    }

    res_test <- switch(which_test,
                       asymptotic = vc_test_asym(y = y_lcpm[genesets, ], x = x, indiv = indiv, phi = phi,
                                                 w = w[genesets, ], Sigma_xi = Sigma_xi,
                                                 genewise_pvals = FALSE, homogen_traj = homogen_traj,
                                                 na.rm = na.rm_tcgsaseq),
                       permutation = vc_test_perm(y = y_lcpm[genesets, ], x = x, indiv = indiv, phi = phi,
                                                  w = w[genesets, ], Sigma_xi = Sigma_xi, n_perm = n_perm,
                                                  genewise_pvals = FALSE, homogen_traj = homogen_traj,
                                                  na.rm = na.rm_tcgsaseq)
    )
    pvals <- data.frame("rawPval" = res_test$set_pval, "adjPval" = NA)
    padjust_methods <- NA

  }

  return(list("which_test" = which_test, "preprocessed" = preprocessed, "n_perm" = n_perm,
              "genesets" = genesets, "pvals" = pvals))

}
