### NOT CURRENTLY IN USE


## A version of limma's voom(), tweaked for single-cell analysis

#' A version of voom() from limma tweaked for single-cell expression data
#'
#' @param counts a numeric matrix containing raw counts, or an
#' ExpressionSet containing raw counts, or a DGEList object.
#' @param design design matrix with rows corresponding to samples and
#' columns to coefficients to be estimated. Defaults to the unit
#' vector meaning that samples are treated as replicates.
#' @param lib.size numeric vector containing total library sizes for
#' each sample. If \code{NULL} and counts is a DGEList then, the
#' normalized library sizes are taken from counts.  Otherwise library
#' sizes are calculated from the columnwise counts totals.
#' @param normalize.method numeric vector containing total library
#' sizes for each sample. If NULL and counts is a DGEList then, the
#' normalized library sizes are taken from counts.  Otherwise library
#' sizes are calculated from the columnwise counts totals.
#' @param plot logical, should a plot of the mean-variance trend be displayed?
#' @param span width of the lowess smoothing window as a proportion.
#' @param prior.count offset, on the scale of raw counts, to be added
#' to raw counts before converting to log2(counts-per-million);
#' default is 0.5.
#' @details This function is intended to process RNA-Seq or ChIP-Seq
#' data prior to linear modelling in limma, written by Charity Law and
#' Gordon Smyth. voom is an acronym for mean-variance modelling at the
#' observational level. The key concern is to estimate the
#' mean-variance relationship in the data, then use this to compute
#' appropriate weights for each observation. Count data almost show
#' non-trivial mean-variance relationships. Raw counts show increasing
#' variance with increasing count size, while log-counts typically
#' show a decreasing mean-variance trend. This function estimates the
#' mean-variance trend for log-counts, then assigns a weight to each
#' observation based on its predicted variance. The weights are then
#' used in the linear modelling process to adjust for
#' heteroscedasticity.  In an experiment, a count value is observed
#' for each tag in each sample. A tag-wise mean-variance trend is
#' computed using lowess. The tag-wise mean is the mean log2 count
#' with a (default) offset of 0.5, across samples for a given tag. The
#' tag-wise variance is the quarter-root-variance of normalized2 log2
#' counts per million values with a (default) offset of 0.5, across
#' samples for a given tag. Tags with zero counts across all samples
#' are not included in the lowess fit. Optional normalization is
#' performed using normalizeBetweenArrays. Using fitted values of log2
#' counts from a linear model fit by lmFit, variances from the
#' mean-variance trend were interpolated for each observation. This
#' was carried out by approxfun. Inverse variance weights can be used
#' to correct for mean-variance trend in the count data.
#' 
#' @return An EList object.
#' @export
#'
# voomsc <- function(counts, design=NULL, lib.size=NULL, normalize.method = "none", plot = FALSE,
#                    span = 0.5, prior.count = 0.5, ...)
#     ## A version of voom() from limma tweaked for single-cell expression data
#     ## Davis McCarthy
#     ## Created 22 October 2014. Last modified 22 October 2014.
#     ## Uses voom() code from:
#     ## Gordon Smyth and Charity Law
#     ## Created 22 June 2011.  Last modified 5 June 2013.
#     ## Linear modelling of count data with mean-variance modelling at the observational level.
#     ## Creates an EList object for entry to lmFit() etc in the limma pipeline.
# {
#     out <- list()
# 
#                                         #	Check counts
#     if(is(counts,"DGEList")) {
#         out$genes <- counts$genes
#         out$targets <- counts$samples
#         if(is.null(design) && diff(range(as.numeric(counts$sample$group)))>0) design <- model.matrix(~group,data=counts$samples)
#         if(is.null(lib.size)) lib.size <- with(counts$samples,lib.size*norm.factors)
#         counts <- counts$counts
#     } else {
#         isExpressionSet <- suppressPackageStartupMessages(is(counts,"ExpressionSet"))
#         if(isExpressionSet) {
#             if(length(fData(counts))) out$genes <- fData(counts)
#             if(length(pData(counts))) out$targets <- pData(counts)
#             counts <- exprs(counts)
#         } else {
#             counts <- as.matrix(counts)
#         }
#     }
# 
#                                         #	Check design
#     if(is.null(design)) {
#         design <- matrix(1,ncol(counts),1)
#         rownames(design) <- colnames(counts)
#         colnames(design) <- "GrandMean"
#     }
# 
#                                         #	Check lib.size
#     if(is.null(lib.size)) lib.size <- colSums(counts)
# 
#                                         #	Fit linear model to log2-counts-per-million
#     y <- t(log2(t(counts + prior.count) / (lib.size+1) * 1e6))
#     y <- normalizeBetweenArrays(y, method = normalize.method)
#     fit <- lmFit(y, design, ...)
#     if(is.null(fit$Amean)) fit$Amean <- rowMeans(y, na.rm=TRUE)
# 
# #	Fit lowess trend to sqrt-standard-deviations by log-count-size
#     sx <- fit$Amean+mean(log2(lib.size + 1)) - log2(1e6)
#     sy <- sqrt(fit$sigma)
#     allzero <- rowSums(counts) == 0
#     if(any(allzero)) {
#         sx <- sx[!allzero]
#         sy <- sy[!allzero]
#     }
#     l <- lowess(sx , sy, f = span)
#     if(plot) {
#         plot(sx, sy, xlab = paste("log2( count size +", prior.count), ylab="Sqrt( standard deviation )", pch = 16, cex = 0.25)
#         title("voom: Mean-variance trend")
#         lines(l,col="red")
#     }
# 
# #	Make interpolating rule
# #	Special treatment of zero counts is now removed;
# #	instead zero counts get same variance as smallest gene average.
# #	l$x <- c(0.5^0.25, l$x)
# #	l$x <- c(log2(0.5), l$x)
# #	var0 <- var(log2(0.5*1e6/(lib.size+0.5)))^0.25
# #	var0 <- max(var0,1e-6)
# #	l$y <- c(var0, l$y)
#     f <- approxfun(l, rule=2)
# 
# #	Find individual quarter-root fitted counts
#     if( is.null(fit$rank) ) {
#         warning("Rank of lmFit is NULL")
#         fit$rank <- ncol(design)
#     }
#     if(fit$rank < ncol(design)) {
#         j <- fit$pivot[1:fit$rank]
#         fitted.values <- fit$coef[,j,drop=FALSE] %*% t(fit$design[,j,drop=FALSE])
#     } else {
#         fitted.values <- fit$coef %*% t(fit$design)
#     }
#     fitted.cpm <- 2^fitted.values
#     fitted.count <- 1e-6 * t(t(fitted.cpm)*(lib.size+1))
#     fitted.logcount <- log2(fitted.count)
# 
# #	Apply trend to individual observations
#     w <- 1/f(fitted.logcount)^4
#     dim(w) <- dim(fitted.logcount)
# 
# #	Output
#     out$E <- y
#     out$weights <- w
#     out$design <- design
#     if(is.null(out$targets))
#         out$targets <- data.frame(lib.size=lib.size)
#     else
#         out$targets$lib.size <- lib.size
#     new("EList", out)
# }
# 
