## Suite of plotting functions

################################################################################
### Generic plot function for SCESet

#' Plot an overview of expression for each cell
#'
#' Plot the relative proportion of the library accounted for by the most highly 
#' expressed features for each cell for an \code{SCESet} dataset.
#' 
#' @param x an \code{SCESet} object
#' @param block1 character string defining the column of \code{pData(object)} to
#' be used as a factor by which to separate the cells into blocks (separate 
#' panels) in the plot. Default is \code{NULL}, in which case there is no 
#' blocking.
#' @param block2 character string defining the column of \code{pData(object)} to
#' be used as a factor by which to separate the cells into blocks (separate 
#' panels) in the plot. Default is \code{NULL}, in which case there is no 
#' blocking.
#' @param colour_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to colour the points in the plot.
#' @param nfeatures numeric scalar indicating the number of features to include 
#' in the plot. 
#' @param use_as_exprs character string indicating which values should be used
#' as the expression values for this plot. Valid arguments are \code{"tpm"} 
#' (default; transcripts per million), \code{"fpkm"} (FPKM values), 
#' \code{"counts"} (counts for each feature) or \code{"exprs"} (whatever is in 
#' the \code{'exprs'} slot of the \code{SCESet} object; if already on the log2
#' scale, as indicated by the \code{logged} slot of the object, then exprs 
#' values are set to the power of 2 (so they are back on the raw scale they were
#'  on) before making the plot).
#' @param linewidth numeric scalar giving the "size" parameter (in ggplot2 
#' parlance) for the lines plotted. Default is 1.5.
#' @param y optional argument for generic \code{plot} functions, not used for 
#' plotting an \code{SCESet} object
#' @param ... arguments passed to \code{plotSCESet}
#' @param ncol number of columns to use for \code{facet_wrap} if only one block is 
#' defined.
#' @param theme_size numeric scalar giving font size to use for the plotting 
#' theme
#' 
#' @details Plots produced by this function are intended to provide an overview
#' of large-scale differences between cells. For each cell, the features are 
#' ordered from most-expressed to least-expressed and the cumulative proportion 
#' of the total expression for the cell is computed across the top 
#' \code{nfeatures} features. These plots can flag cells with a very high 
#' proportion of the library coming from a small number of features; such cells 
#' are likely to be problematic for analyses. Using the colour and blocking 
#' arguments can flag overall differences in cells under different experimental
#' conditions or affected by different batch and other variables.
#' 
#' @name plot
#' @aliases plot plot,SCESet-method plot,SCESet,ANY-method
#' @export
#' @method plot
#'
#' @examples
#' ## Set up an example SCESet
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' 
#' plot(example_sceset, use_as_exprs = "exprs")
#' plot(example_sceset, use_as_exprs = "exprs", colour_by = "Cell_Cycle")
#' plot(example_sceset, use_as_exprs = "exprs", block1 = "Treatment", 
#' colour_by = "Cell_Cycle")
#' plot(example_sceset, use_as_exprs = "exprs", block1 = "Treatment", 
#' block2 = "Mutation_Status", colour_by = "Cell_Cycle")
#' # What happens if chosen expression values are not available?
#' plot(example_sceset, block1 = "Treatment", colour_by = "Cell_Cycle") 
#' 
setMethod("plot", signature("SCESet"),
          function(x, ...) {
              plotSCESet(x, ...)
          })

#' @rdname plot
#' @aliases plot
#' @export
plotSCESet <- function(x, block1 = NULL, block2 = NULL, colour_by = NULL, 
                        nfeatures = 500, use_as_exprs = "tpm", ncol = 3, 
                       linewidth = 1.5, theme_size = 10) {
    object <- x
    if ( !is(object, "SCESet") )
        stop("Object must be an SCESet")
    if ( !is.null(block1) ) {
        if ( !(block1 %in% colnames(pData(object))) )
            stop("The block1 argument must either be NULL or a column of pData(object).")        
    } 
    if ( !is.null(block2) ) {
        if ( !(block2 %in% colnames(pData(object))) )
            stop("The block2 argument must either be NULL or a column of pData(object).")        
    }  
    if ( !is.null(colour_by) ) {
        if ( !(colour_by %in% colnames(pData(object))) )
            stop("The colour_by argument must either be NULL or a column of pData(object).")        
    }  
    
    ## Define an expression matrix depending on which values we're using
    use_as_exprs <- match.arg(use_as_exprs, c("exprs", "tpm", "fpkm", "counts"))
    exprs_mat <- switch(use_as_exprs,
                        exprs = exprs(object),
                        tpm = tpm(object),
                        fpkm = fpkm(object),
                        counts = counts(object))
    if ( is.null(exprs_mat) ) {
        warning(paste0("The object does not contain ", use_as_exprs, " expression values. Using exprs(object) values instead."))
        exprs_mat <- exprs(object)
        use_as_exprs <- "exprs"
    }
    if ( use_as_exprs == "exprs" && object@logged )
        exprs_mat <- 2 ^ (exprs_mat) - object@logExprsOffset
    
    ## Use plyr to get the sequencing real estate accounted for by features
    nfeatures_total <- nrow(exprs_mat)
    seq_real_estate <- t(plyr::aaply(exprs_mat, 2, .fun = function(x) {
        cumsum(sort(x, decreasing = TRUE))
    }))
    rownames(seq_real_estate) <- 1:nfeatures_total    
    nfeatures_to_plot <- nfeatures
    seq_real_estate_long <- reshape2::melt(seq_real_estate[1:nfeatures_to_plot, ], 
                                           value.name = "exprs")
    
    ## Get the proportion of the library accounted for by the top features
    prop_library <- reshape2::melt(t(t(seq_real_estate[1:nfeatures_to_plot, ]) /
                                         colSums(exprs_mat)), 
                                   value.name = "prop_library")
    colnames(seq_real_estate_long) <- c("Feature", "Cell", "exprs")
    seq_real_estate_long$Proportion_Library <- prop_library$prop_library
    
    ## Add block and colour_by information if provided
    if ( !is.null(block1) )
        seq_real_estate_long <- dplyr::mutate(
            seq_real_estate_long, block1 = as.factor(rep(object[[block1]], 
                                                       each = nfeatures_to_plot)))
    if ( !is.null(block2) )                                  
        seq_real_estate_long <- dplyr::mutate(
            seq_real_estate_long, block2 = as.factor(rep(object[[block2]], 
                                                       each = nfeatures_to_plot)))
    if ( !is.null(colour_by) )                                  
        seq_real_estate_long <- dplyr::mutate(
            seq_real_estate_long, colour_by = rep(object[[colour_by]], 
                                             each = nfeatures_to_plot))
    ## Set up plot
    if ( is.null(colour_by) ) {
        plot_out <- ggplot(seq_real_estate_long, 
                           aes_string(x = "Feature", y = "Proportion_Library", 
                                      group = "Cell")) +
            geom_line(linetype = "solid", alpha = 0.3, size = linewidth)
    } else {
        plot_out <- ggplot(seq_real_estate_long, 
                           aes_string(x = "Feature", y = "Proportion_Library", 
                                      group = "Cell", colour = "colour_by")) +
            geom_line(linetype = "solid", alpha = 0.3, size = linewidth)
    }
    ## Deal with blocks for grid
    if ( !(is.null(block1) | is.null(block2)) )
        plot_out <- plot_out + facet_grid(block2 ~ block1)
    else {
        if ( !is.null(block1) & is.null(block2) ) {
            plot_out <- plot_out + 
                facet_wrap(~block1, ncol = ncol)
        }
        if ( is.null(block1) & !is.null(block2) ) {
            plot_out <- plot_out + 
                facet_wrap(~block2, ncol = ncol)
        }
    }
    ## Add extra plot theme and details
    plot_out <- plot_out + ggthemes::scale_colour_tableau(name = colour_by) + 
        xlab("Number of features") + ylab("Cumulative proportion of library")
    
    if ( library(cowplot, logical.return = TRUE) )
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    else
        plot_out <- plot_out + theme_bw(theme_size)
    ## Return plot
    plot_out
}


################################################################################

#' Plot PCA for an SCESet object
#'
#' Produce a principal components analysis (PCA) plot of two or more principal 
#' components for an \code{SCESet} dataset.
#' 
#' @param object an \code{SCESet} object
#' @param ntop numeric scalar indicating the number of most variable features to 
#' use for the PCA. Default is \code{500}, but any \code{ntop} argument is 
#' overrided if the \code{feature_set} argument is non-NULL.
#' @param ncomponents numeric scalar indicating the number of principal 
#' components to plot, starting from the first principal component. Default is 
#' 2. If \code{ncomponents} is 2, then a scatterplot of PC2 vs PC1 is produced.
#' If \code{ncomponents} is greater than 2, a pairs plots for the top components
#' is produced.
#' @param colour_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to colour the points in the plot.
#' @param shape_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to define the shape of the points in the plot.
#' @param size_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to define the size of points in the plot.
#' @param feature_set character, numeric or logical vector indicating a set of
#' features to use for the PCA. If character, entries must all be in 
#' \code{featureNames(object)}. If numeric, values are taken to be indices for 
#' features. If logical, vector is used to index features and should have length
#' equal to \code{nrow(object)}.
#' @param return_SCESet logical, should the function return an \code{SCESet} 
#' object with principal component values for cells in the 
#' \code{reducedDimension} slot. Default is \code{FALSE}, in which case a 
#' \code{ggplot} object is returned.
#' @param scale_features logical, should the expression values be standardised
#' so that each feature has unit variance? Default is \code{TRUE}.
#' @param draw_plot logical, should the plot be drawn on the current graphics 
#' device? Only used if \code{return_SCESet} is \code{TRUE}, otherwise the plot 
#' is always produced.
#' @param theme_size numeric scalar giving default font size for plotting theme 
#' (default is 10).
#' 
#' @details The function \code{\link{prcomp}} is used internally to do the PCA. 
#' The function checks whether the \code{object} has standardised 
#' expression values (by looking at \code{stand_exprs(object)}). If yes, the 
#' existing standardised expression values are used for the PCA. If not, then 
#' standardised expression values are computed using \code{\link{scale}} (with
#' feature-wise unit variances or not according to the \code{scale_features} 
#' argument), added to the object and PCA is done using these new standardised
#' expression values.
#'
#' @name plotPCA
#' @aliases plotPCA plotPCA,SCESet-method
#' @export
#'
#' @examples
#' ## Set up an example SCESet
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' 
#' ## Examples plotting PC1 and PC2
#' plotPCA(example_sceset)
#' plotPCA(example_sceset, colour_by = "Cell_Cycle")
#' plotPCA(example_sceset, colour_by = "Cell_Cycle", shape_by = "Treatment")
#' plotPCA(example_sceset, colour_by = "Cell_Cycle", shape_by = "Treatment", 
#' size_by = "Mutation_Status")
#' plotPCA(example_sceset, shape_by = "Treatment", size_by = "Mutation_Status")
#' plotPCA(example_sceset, feature_set = 1:100, colour_by = "Treatment", 
#' shape_by = "Mutation_Status")
#' 
#' plotPCA(example_sceset, shape_by = "Treatment", return_SCESet = TRUE)
#' 
#' ## Examples plotting more than 2 PCs
#' plotPCA(example_sceset, ncomponents = 8)
#' plotPCA(example_sceset, ncomponents = 4, colour_by = "Treatment", 
#' shape_by = "Mutation_Status")
#' 
plotPCASCESet <- function(object, ntop=500, ncomponents=2, colour_by=NULL, 
                          shape_by=NULL, size_by=NULL, feature_set=NULL, 
                          return_SCESet=FALSE, scale_features=TRUE, 
                          draw_plot = TRUE, theme_size = 10) {
    ## Check arguments are valid
    if ( !is.null(colour_by) ) {
        if ( !(colour_by %in% varLabels(object)) )
            stop("the argument 'colour_by' should specify a column of pData(object) [see varLabels(object)]")
    } else {
        if ( "is_cell_control" %in% varLabels(object) )
            colour_by <- "is_cell_control"
    }
    if ( !is.null(shape_by) ) {
        if ( !(shape_by %in% varLabels(object)) )
            stop("the argument 'shape_by' should specify a column of pData(object) [see varLabels(object)]")
        if ( nlevels(as.factor(pData(object)[[shape_by]])) > 10 )
            stop("when coerced to a factor, 'shape_by' should have fewer than 10 levels")   
    } else {
        if ( "is_cell_control" %in% varLabels(object) )
            shape_by <- "is_cell_control"
    }
    if ( !is.null(size_by) ) {
        if ( !(size_by %in% varLabels(object)) )
            stop("the argument 'size_by' should specify a column of pData(object) [see varLabels(object)]")
    }
    if ( !is.null(feature_set) && typeof(feature_set) == "character" ) {
        if ( !(all(feature_set %in% featureNames(object))) )
            stop("when the argument 'feature_set' is of type character, all features must be in featureNames(object)")
    }
    
    ## Define features to use: either ntop, or if a set of features is defined, then those
    if ( is.null(feature_set) ) {
        rv <- matrixStats::rowVars(exprs(object))
        feature_set <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    }
#     ## Standardise expression if stand_exprs(object) is null
#     if ( is.null(stand_exprs(object)) ) {
#         stand_exprs(object) <- t(scale(t(exprs(object)), scale = scale_features)) 
#         message("stand_exprs(object) was null, so standardising exprs values for PCA.")
#     } else
#         message("stand_exprs(object) was not null, so using these values for PCA.")

    ## Standardise expression if stand_exprs(object) is null
    exprs_to_plot <- t(scale(t(exprs(object)), scale = scale_features)) 
    
    ## Drop any features with zero variance
    exprs_to_plot <- exprs_to_plot[feature_set,]
    keep_feature <- (matrixStats::rowVars(exprs_to_plot) > 0.001)
    keep_feature[is.na(keep_feature)] <- FALSE
    exprs_to_plot <- exprs_to_plot[keep_feature, ]
    
    ## Compute PCA
    pca <- prcomp(t(exprs_to_plot))
    percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
    
    ## Define data.frame for plotting
    df_to_plot <- data.frame(pca$x[, 1:ncomponents], 
                             row.names = sampleNames(object))
    if ( !is.null(colour_by) ) {
        if ( nlevels(as.factor(pData(object)[[colour_by]])) > 10 )
            df_to_plot <- data.frame(
                df_to_plot, colour_by = pData(object)[[colour_by]])
        else
            df_to_plot <- data.frame(
                df_to_plot, colour_by = as.factor(pData(object)[[colour_by]]))
    }
    if ( !is.null(shape_by) )
        df_to_plot <- data.frame(df_to_plot, 
                                 shape_by = as.factor(pData(object)[[shape_by]]))
    if ( !is.null(size_by) )
        df_to_plot <- data.frame(df_to_plot, size_by = pData(object)[[size_by]])
    
    plot_out <- plotReducedDim.default(df_to_plot, ncomponents, colour_by, shape_by, 
                                       size_by, percentVar)
    
    ## Define plotting theme
    if ( library(cowplot, logical.return = TRUE) )
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    else
        plot_out <- plot_out + theme_bw(theme_size)
    
    ## Plot PCA and return appropriate object
    if (return_SCESet) {
        df_out <- pca$x[, 1:ncomponents]
        rownames(df_out) <- sampleNames(object)
        attr(df_out, "percentVar") <- percentVar[1:ncomponents]
        reducedDimension(object) <- df_out
        if ( draw_plot )
            print(plot_out)
        return(object)
    } else {
        ## Return PCA plot
        return(plot_out)
    }
}
    

#' @rdname plotPCA
#' @aliases plotPCA
#' @export
setMethod("plotPCA", signature("SCESet"),
          function(object, ntop=500, ncomponents=2, colour_by=NULL, shape_by=NULL, 
                   size_by=NULL, feature_set=NULL, return_SCESet=FALSE, 
                   scale_features=TRUE, draw_plot = TRUE, theme_size = 10) {
              plotPCASCESet(object, ntop, ncomponents, colour_by, shape_by, size_by, 
                             feature_set, return_SCESet, scale_features, draw_plot, theme_size)
          })


.makePairs <- function(data_matrix) {
    ## with thanks to Gaston Sanchez, who posted this code online
    ## https://gastonsanchez.wordpress.com/2012/08/27/scatterplot-matrices-with-ggplot/
    if ( is.null(names(data_matrix)) )
        names(data_matrix) <- paste0("row", 1:nrow(data_matrix))
    exp_grid <- expand.grid(x = 1:ncol(data_matrix), y = 1:ncol(data_matrix))
    exp_grid <- subset(exp_grid, x != y)
    all_panels <- do.call("rbind", lapply(1:nrow(exp_grid), function(i) {
        xcol <- exp_grid[i, "x"]
        ycol <- exp_grid[i, "y"]
        data.frame(xvar = names(data_matrix)[ycol], yvar = names(data_matrix)[xcol], 
                   x = data_matrix[, xcol], y = data_matrix[, ycol], data_matrix)
    }))
    all_panels$xvar <- factor(all_panels$xvar, levels = names(data_matrix))
    all_panels$yvar <- factor(all_panels$yvar, levels = names(data_matrix))
    densities <- do.call("rbind", lapply(1:ncol(data_matrix), function(i) {
        data.frame(xvar = names(data_matrix)[i], yvar = names(data_matrix)[i], 
                   x = data_matrix[, i])
    }))
    list(all = all_panels, densities = densities)
}


################################################################################
### plotTSNE

#' Plot t-SNE for an SCESet object
#'
#' Produce a t-distributed stochastic neighbour embedding plot of two components
#'  for an \code{SCESet} dataset.
#' 
#' @param object an \code{SCESet} object
#' @param ntop numeric scalar indicating the number of most variable features to 
#' use for the PCA. Default is \code{500}, but any \code{ntop} argument is 
#' overrided if the \code{feature_set} argument is non-NULL.
#' @param ncomponents numeric scalar indicating the number of principal 
#' components to plot, starting from the first principal component. Default is 
#' 2. If \code{ncomponents} is 2, then a scatterplot of PC2 vs PC1 is produced.
#' If \code{ncomponents} is greater than 2, a pairs plots for the top components
#' is produced. NB: computing more than two components for t-SNE can become very
#' time consuming.
#' @param colour_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to colour the points in the plot.
#' @param shape_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to define the shape of the points in the plot.
#' @param size_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to define the size of points in the plot.
#' @param feature_set character, numeric or logical vector indicating a set of
#' features to use for the PCA. If character, entries must all be in 
#' \code{featureNames(object)}. If numeric, values are taken to be indices for 
#' features. If logical, vector is used to index features and should have length
#' equal to \code{nrow(object)}.
#' @param return_SCESet logical, should the function return an \code{SCESet} 
#' object with principal component values for cells in the 
#' \code{reducedDimension} slot. Default is \code{FALSE}, in which case a 
#' \code{ggplot} object is returned.
#' @param scale_features logical, should the expression values be standardised
#' so that each feature has unit variance? Default is \code{TRUE}.
#' @param draw_plot logical, should the plot be drawn on the current graphics 
#' device? Only used if \code{return_SCESet} is \code{TRUE}, otherwise the plot 
#' is always produced.
#' @param theme_size numeric scalar giving default font size for plotting theme 
#' (default is 10).
#' @param rand_seed (optional) numeric scalar that can be passed to 
#' \code{set.seed} to make plots reproducible.
#' @param perplexity numeric scalar value defining the "perplexity parameter" 
#' for the t-SNE plot. Passed to \code{\link[Rtsne]{Rtsne}} - see documentation
#' for that package for more details.
#' @param ... further arguments passed to \code{\link[Rtsne]{Rtsne}}
#' 
#' @details The function \code{\link[Rtsne]{Rtsne}} is used internally to 
#' compute the t-SNE. The function checks whether the \code{object} has 
#' standardised expression values (by looking at \code{stand_exprs(object)}). If
#' yes, the existing standardised expression values are used for the PCA. If 
#' not, then standardised expression values are computed using 
#' \code{\link{scale}} (with feature-wise unit variances or not according to the
#' \code{scale_features} argument), added to the object and PCA is done using 
#' these new standardised expression values.
#'
#' @return If \code{return_SCESet} is \code{TRUE}, then the function returns an
#' \code{SCESet} object, otherwise it returns a \code{ggplot} object.
#' @name plotTSNE
#' @aliases plotTSNE plotTSNE,SCESet-method
#' @import Rtsne
#' @export
#' @seealso 
#' \code{\link[Rtsne]{Rtsne}}
#' @references 
#' L.J.P. van der Maaten. Barnes-Hut-SNE. In Proceedings of the International 
#' Conference on Learning Representations, 2013.
#' 
#' @examples
#' ## Set up an example SCESet
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' 
#' ## Examples plotting PC1 and PC2
#' plotTSNE(example_sceset, perplexity = 10)
#' plotTSNE(example_sceset, colour_by = "Cell_Cycle", perplexity = 10)
#' plotTSNE(example_sceset, colour_by = "Cell_Cycle", shape_by = "Treatment",
#' perplexity = 10)
#' plotTSNE(example_sceset, colour_by = "Cell_Cycle", shape_by = "Treatment", 
#' size_by = "Mutation_Status", perplexity = 10)
#' plotTSNE(example_sceset, shape_by = "Treatment", size_by = "Mutation_Status",
#' perplexity = 5)
#' plotTSNE(example_sceset, feature_set = 1:100, colour_by = "Treatment", 
#' shape_by = "Mutation_Status", perplexity = 5)
#' 
#' plotTSNE(example_sceset, shape_by = "Treatment", return_SCESet = TRUE,
#' perplexity = 10)
#' 
#' 
setMethod("plotTSNE", signature("SCESet"),
          function(object, ntop = 500, ncomponents = 2, colour_by = NULL, 
                   shape_by = NULL, size_by = NULL, feature_set = NULL, 
                   return_SCESet = FALSE, scale_features = TRUE, 
                   draw_plot = TRUE, theme_size = 10, rand_seed = NULL, 
                   perplexity = floor(ncol(object) / 5), ...) {
              ##
              if ( !library(Rtsne, logical.return = TRUE) )
                  stop("This function requires the 'Rtsne' package. 
                       Try: install.packages('Rtsne').")
              ## Check arguments are valid
              if ( !is.null(colour_by) ) {
                  if ( !(colour_by %in% varLabels(object)) )
                      stop("the argument 'colour_by' should specify a column of 
                           pData(object) [see varLabels(object)]")
              } else {
                  if ( "is_cell_control" %in% varLabels(object) )
                      colour_by <- "is_cell_control"
              }
              if ( !is.null(shape_by) ) {
                  if ( !(shape_by %in% varLabels(object)) )
                      stop("the argument 'shape_by' should specify a column of 
                           pData(object) [see varLabels(object)]")
                  if ( nlevels(as.factor(pData(object)[[shape_by]])) > 10 )
                      stop("when coerced to a factor, 'shape_by' should have 
                           fewer than 10 levels")   
              } else {
                  if ( "is_cell_control" %in% varLabels(object) )
                      shape_by <- "is_cell_control"
              }
              if ( !is.null(size_by) ) {
                  if ( !(size_by %in% varLabels(object)) )
                      stop("the argument 'size_by' should specify a column of 
                           pData(object) [see varLabels(object)]")
              }
              if ( !is.null(feature_set) && typeof(feature_set) == "character" ) {
                  if ( !(all(feature_set %in% featureNames(object))) )
                      stop("when the argument 'feature_set' is of type character,
                           all features must be in featureNames(object)")
              }
              
              ## Define features to use: either ntop, or if a set of features is
              ## defined, then those
              if ( is.null(feature_set) ) {
                  rv <- matrixStats::rowVars(exprs(object))
                  feature_set <- 
                      order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                               length(rv)))]
              }
              ## Standardise expression if stand_exprs(object) is null
              exprs_to_plot <- t(scale(t(exprs(object)), scale = scale_features)) 
              
              ## Drop any features with zero variance
              exprs_to_plot <- exprs_to_plot[feature_set,]
              keep_feature <- (matrixStats::rowVars(exprs_to_plot) > 0.001)
              keep_feature[is.na(keep_feature)] <- FALSE
              exprs_to_plot <- exprs_to_plot[keep_feature, ]
              
              ## Compute t-SNE
              if ( !is.null(rand_seed) )
                  set.seed(rand_seed)
              tsne_out <- Rtsne::Rtsne(t(exprs_to_plot),
                                       initial_dims = max(50, ncol(object)),
                                       perplexity = perplexity, ...)
              
              
              ## Define data.frame for plotting
              df_to_plot <- data.frame(tsne_out$Y[, 1:ncomponents], 
                                       row.names = sampleNames(object))
              if ( !is.null(colour_by) ) {
                  if ( nlevels(as.factor(pData(object)[[colour_by]])) > 10 )
                      df_to_plot <- data.frame(
                          df_to_plot, colour_by = pData(object)[[colour_by]])
                  else
                      df_to_plot <- data.frame(
                          df_to_plot, 
                          colour_by = as.factor(pData(object)[[colour_by]]))
              }
              if ( !is.null(shape_by) )
                  df_to_plot <- data.frame(
                      df_to_plot, 
                      shape_by = as.factor(pData(object)[[shape_by]]))
              if ( !is.null(size_by) )
                  df_to_plot <- data.frame(df_to_plot, 
                                           size_by = pData(object)[[size_by]])
              
              plot_out <- plotReducedDim.default(df_to_plot, ncomponents, 
                                                 colour_by, shape_by, size_by)
              
              ## Define plotting theme
              if ( library(cowplot, logical.return = TRUE) )
                  plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
              else
                  plot_out <- plot_out + theme_bw(theme_size)
              
              ## Plot PCA and return appropriate object
              if (return_SCESet) {
                  df_out <- tsne_out$Y[, 1:ncomponents]
                  rownames(df_out) <- sampleNames(object)
                  reducedDimension(object) <- df_out
                  if ( draw_plot )
                      print(plot_out)
                  return(object)
              } else {
                  ## Return PCA plot
                  return(plot_out)
              }
          })


################################################################################
### plotReducedDim

#' Plot reduced dimension representation of cells
#' 
#' @param object an \code{SCESet} object with a non-NULL \code{reducedDimension}
#' slot.
#' @param df_to_plot data.frame containing a reduced dimension represenation of
#' cells and optional metadata for the plot.
#' @param ncomponents numeric scalar indicating the number of principal 
#' components to plot, starting from the first principal component. Default is 
#' 2. If \code{ncomponents} is 2, then a scatterplot of Dimension 2 vs Dimension 
#' 1 is produced. If \code{ncomponents} is greater than 2, a pairs plots for the
#'  top dimensions is produced.
#' @param colour_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to colour the points in the plot.
#' @param shape_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to define the shape of the points in the plot.
#' @param size_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to define the size of points in the plot.
#' @param percentVar numeric vector giving the proportion of variance in 
#' expression explained by each reduced dimension. Only expected to be used 
#' internally in the \code{\link[scater]{plotPCA}} function.
#' @param theme_size numeric scalar giving default font size for plotting theme 
#' (default is 10).
#' @param ... optional arguments (from those listed above) passed to 
#' \code{plotReducedDim.SCESet} or \code{plotReducedDim.default}
#' 
#' @details The function \code{plotReducedDim.default} assumes that the first 
#' \code{ncomponents} columns of \code{df_to_plot} contain the reduced dimension
#'  components to plot, and that any subsequent columns define factors for 
#'  \code{colour_by}, \code{shape_by} and \code{size_by} in the plot.
#' 
#' @name plotReducedDim
#' @aliases plotReducedDim plotReducedDim,SCESet-method plotReducedDim,data.frame-method 
#' @export
#' 
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' reducedDimension(example_sceset) <- prcomp(t(exprs(example_sceset)), scale. = TRUE)$x
#' plotReducedDim(example_sceset)
#' plotReducedDim(example_sceset, colour_by="Cell_Cycle")
#' plotReducedDim(example_sceset, colour_by="Cell_Cycle", shape_by="Treatment")
#' plotReducedDim(example_sceset, colour_by="Cell_Cycle", size_by="Treatment")
#' plotReducedDim(example_sceset, ncomponents=5)
#' plotReducedDim(example_sceset, ncomponents=5, colour_by="Cell_Cycle", shape_by="Treatment")
#' 
#' 
plotReducedDim.default <- function(df_to_plot, ncomponents=2, colour_by=NULL, 
                           shape_by=NULL, size_by=NULL, percentVar=NULL,
                           theme_size = 10) {
    ## Define plot
    if ( ncomponents > 2 ) {
        ## expanding numeric columns for pairs plot
        df_to_expand <- df_to_plot[, 1:ncomponents]
        if ( is.null(percentVar) ) {
            colnames(df_to_expand) <- colnames(df_to_plot)[1:ncomponents]
        } else {
            colnames(df_to_expand) <- paste0(colnames(df_to_plot)[1:ncomponents], ": ", 
                                             round(percentVar[1:ncomponents] * 100), 
                                             "% variance")
        }
        gg1 <- .makePairs(df_to_expand)
        ## new data frame 
        df_to_plot_big <- data.frame(gg1$all, df_to_plot[, -c(1:ncomponents)])
        colnames(df_to_plot_big)[-c(1:4)] <- colnames(df_to_plot)
        ## pairs plot
        plot_out <- ggplot(df_to_plot_big, aes_string(x = "x", y = "y")) + 
            facet_grid(xvar ~ yvar, scales = "free") + 
            stat_density(aes_string(x = "x", 
                                    y = "(..scaled.. * diff(range(x)) + min(x))"),
                         data = gg1$densities, position = "identity", 
                         colour = "grey20", geom = "line") +
            xlab("") + 
            ylab("") +
            theme_bw(theme_size)
    } else {
        comps <- colnames(df_to_plot)[1:2]
        if ( is.null(percentVar) ) {
            x_lab <- "Dimension 1"
            y_lab <- "Dimension 2"
        } else {
            x_lab <- paste0("Component 1: ", round(percentVar[1] * 100), "% variance")
            y_lab <- paste0("Component 2: ", round(percentVar[2] * 100), "% variance")
        }
        plot_out <- ggplot(df_to_plot, aes_string(x = comps[1], y = comps[2])) + 
            xlab(x_lab) + 
            ylab(y_lab) +
            geom_rug(colour = "gray20", alpha = 0.65) +
            theme_bw(theme_size)        
    }
    
    ## Apply colour_by, shape_by and size_by variables if defined
    if ( !is.null(colour_by) & !is.null(shape_by) & !is.null(size_by) ) {
        plot_out <- plot_out + 
            geom_point(aes_string(colour = "colour_by", shape = "shape_by", 
                                  size = "size_by"), alpha = 0.65) +
            guides(size = guide_legend(title = size_by), 
                   shape = guide_legend(title = shape_by))
        if ( is.numeric(df_to_plot$colour_by) ) {
            plot_out <- plot_out + scale_colour_gradient(
                name = colour_by, low = "gold", high = "darkred", space = "Lab")
        } else {
            plot_out <- plot_out + 
                ggthemes::scale_colour_tableau(name = colour_by)
        }
    } else {
        if  ( sum(is.null(colour_by) + is.null(shape_by) + is.null(size_by)) == 1 ) {
            if ( !is.null(colour_by) & !is.null(shape_by) ) {
                plot_out <- plot_out + 
                    geom_point(aes_string(colour = "colour_by", 
                                          shape = "shape_by"), size = 4, 
                               alpha = 0.65) +
                    guides(shape = guide_legend(title = shape_by))
                if ( is.numeric(df_to_plot$colour_by) ) {
                    plot_out <- plot_out + scale_colour_gradient(
                        name = colour_by, low = "gold", high = "darkred", 
                        space = "Lab")
                } else {
                    plot_out <- plot_out + 
                        ggthemes::scale_colour_tableau(name = colour_by)
                }
            }
            if ( !is.null(colour_by) & !is.null(size_by) ) {
                plot_out <- plot_out + 
                    geom_point(aes_string(fill = "colour_by", size = "size_by"), 
                               shape = 21, colour = "gray70", alpha = 0.65) +
                    guides(size = guide_legend(title = size_by))
                if ( is.numeric(df_to_plot$colour_by) ) {
                    plot_out <- plot_out + scale_colour_gradient(
                        name = colour_by, low = "gold", high = "darkred", 
                        space = "Lab")
                } else {
                    plot_out <- plot_out + 
                        ggthemes::scale_colour_tableau(name = colour_by)
                }
            }
            if ( !is.null(shape_by) & !is.null(size_by) ) {
                plot_out <- plot_out + 
                    geom_point(aes_string(shape = "shape_by", size = "size_by"), 
                               fill = "gray20", colour = "gray20", alpha = 0.65) +
                    guides(size = guide_legend(title = size_by), 
                           shape = guide_legend(title = shape_by))
            } 
        } else {
            if ( sum(is.null(colour_by) + is.null(shape_by) + is.null(size_by)) == 2 ) {
                if ( !is.null(colour_by) ) {
                    plot_out <- plot_out + 
                        geom_point(aes_string(fill = "colour_by"), size = 4, 
                                   shape = 21, colour = "gray70", alpha = 0.65)
                    if ( is.numeric(df_to_plot$colour_by) ) {
                        plot_out <- plot_out + scale_colour_gradient(
                            name = colour_by, low = "gold", high = "darkred", 
                            space = "Lab")
                    } else {
                        plot_out <- plot_out + 
                            ggthemes::scale_colour_tableau(name = colour_by)
                    }
                }
                if ( !is.null(shape_by) ) {
                    plot_out <- plot_out + 
                        geom_point(aes_string(shape = "shape_by"), size = 4, 
                                   colour = "gray20", alpha = 0.65) +
                        guides(shape = guide_legend(title = shape_by))
                }
                if ( !is.null(size_by) ) {
                    plot_out <- plot_out + 
                        geom_point(aes_string(size = "size_by"), fill = "gray20",
                                   shape = 21, colour = "gray70", alpha = 0.65) +
                        guides(size = guide_legend(title = size_by))
                } 
            } else {
                plot_out <- plot_out + 
                    geom_point(size = 4, fill = "gray20", shape = 21, 
                               colour = "gray70", alpha = 0.65)
            }            
        }
    }

    ## Return plot
    plot_out
}

#' @rdname plotReducedDim
#' @aliases plotReducedDim
#' @export
plotReducedDim.SCESet <- function(object, ncomponents=2, colour_by=NULL, 
                                  shape_by=NULL, size_by=NULL, theme_size = 10) {
    ## Check arguments are valid
    if ( !is.null(colour_by) ) {
        if ( !(colour_by %in% varLabels(object)) )
            stop("the argument 'colour_by' should specify a column of pData(object) [see varLabels(object)]")
    }
    if ( !is.null(shape_by) ) {
        if ( !(shape_by %in% varLabels(object)) )
            stop("the argument 'shape_by' should specify a column of pData(object) [see varLabels(object)]")
        if ( nlevels(as.factor(pData(object)[[shape_by]])) > 10 )
            stop("when coerced to a factor, 'shape_by' should have fewer than 10 levels")   
    }
    if ( !is.null(size_by) ) {
        if ( !(size_by %in% varLabels(object)) )
            stop("the argument 'size_by' should specify a column of pData(object) [see varLabels(object)]")
    }
    
    ## Extract reduced dimension representation of cells
    if ( is.null(reducedDimension(object)) )
        stop("reducedDimension slot of object is NULL. Need non null reducedDimension to plot.")
    red_dim <- redDim(object)
    if ( ncomponents > ncol(red_dim) )
        stop("ncomponents to plot is larger than number of columns of reducedDimension(object)")
    
    ## Define data.frame for plotting
    df_to_plot <- data.frame(red_dim[, 1:ncomponents])
    if ( !is.null(colour_by) ) {
        if ( nlevels(as.factor(pData(object)[[colour_by]])) > 10 )
            df_to_plot <- data.frame(
                df_to_plot, colour_by = pData(object)[[colour_by]])
        else
            df_to_plot <- data.frame(
                df_to_plot, colour_by = as.factor(pData(object)[[colour_by]]))
    }
    if ( !is.null(shape_by) )
        df_to_plot <- data.frame(df_to_plot, 
                                 shape_by = as.factor(pData(object)[[shape_by]]))
    if ( !is.null(size_by) )
        df_to_plot <- data.frame(df_to_plot, size_by = pData(object)[[size_by]])
    
    ## Call default method to make the plot
    plot_out <- plotReducedDim.default(df_to_plot, ncomponents, colour_by, 
                                       shape_by, size_by, percentVar = NULL,
                                       theme_size)
    
    ## Define plotting theme
    if ( library(cowplot, logical.return = TRUE) )
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    else
        plot_out <- plot_out + theme_bw(theme_size)
    plot_out
}

#' @rdname plotReducedDim
#' @aliases plotReducedDIm
#' @export
setMethod("plotReducedDim", signature("SCESet"),
          function(object, ncomponents=2, colour_by=NULL, shape_by=NULL, 
                   size_by=NULL) {
              plotReducedDim.SCESet(object, ncomponents, colour_by, shape_by, 
                                    size_by)
          })

#' @rdname plotReducedDim
#' @aliases plotReducedDIm
#' @export
setMethod("plotReducedDim", signature("data.frame"),
          function(object, ncomponents=2, colour_by=NULL, shape_by=NULL, 
                   size_by=NULL, percentVar=NULL) {
              plotReducedDim.default(object, ncomponents, colour_by, shape_by, 
                                     size_by, percentVar)
          })


################################################################################

#' Plot expression values for a set of features (e.g. genes or transcripts)
#'
#' @param object an SCESet object containing expression values and
#' experimental information. Must have been appropriately prepared. For the 
#' \code{plotExpressionDefault} method, the \code{object} argument is a 
#' data.frame in 'long' format providing expression values for a set of features
#'  to plot, plus metadata used in the \code{aesth} argument, but this is not
#'  meant to be a user-level operation.
#' @param features a character vector of feature names or Boolean
#' vector or numeric vector of indices indicating which features should have 
#' their expression values plotted
#' @param x character string providing a column name of \code{pData(object)} to 
#' plot on the x-axis in the expression plot(s)
#' @param use_as_exprs character string indicating which expression values to plot:
#' either "exprs" for the expression value defined in the \code{exprs} slot or
#' "counts" to plot log2(counts-per-million + 1) using counts from the 
#' \code{counts} slot of \code{object}
#' @param colour_by optional character string supplying name of a column of 
#' \code{pData(object)} which will be used as a variable by which to colour 
#' expression values on the plot.
#' @param shape_by optional character string supplying name of a column of 
#' \code{pData(object)} which will be used as a variable to define the shape of
#' points for expression values on the plot.
#' @param size_by optional character string supplying name of a column of 
#' \code{pData(object)} which will be used as a variable to define the size of
#' points for expression values on the plot.
#' @param ncol number of columns to be used for the panels of the plot
#' @param xlab label for x-axis; if \code{NULL} (default), then \code{x} will be
#' used as the x-axis label
#' @param show_median logical, show the median for each group on the plot
#' @param show_violin logical, show a violin plot for the distribution
#' for each group on the plot
#' @param theme_size numeric scalar giving default font size for plotting theme 
#' (default is 10)
#' @param ... optional arguments (from those listed above) passed to 
#' \code{plotExpressionSCESet} or \code{plotExpressionDefault}
#'
#' @details Plot expression values (default log2(transcripts-per-million +
#' 1), if available) for a set of features.
#' 
#' @name plotExpression
#' @aliases plotExpression plotExpression,SCESet-method plotExpression,data.frame-method
#' @import ggplot2
#' @importFrom Biobase varLabels
#' @importFrom Biobase featureNames
#' @export
#' 
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset <- calculateQCMetrics(example_sceset)
#' plotExpression(example_sceset, 1:6, x="Mutation_Status", use_as_exprs="exprs", 
#' colour_by="Cell_Cycle", show_violin=TRUE, show_median=TRUE)
#' plotExpression(example_sceset, 1:6, x="Mutation_Status", use_as_exprs="counts", 
#' colour_by="Cell_Cycle", show_violin=TRUE, show_median=TRUE)
#' 
plotExpressionSCESet <- function(object, features, x, use_as_exprs = "exprs", 
                                 colour_by = NULL, shape_by = NULL, 
                                 size_by = NULL, ncol = 2, xlab = NULL, 
                                 show_median = TRUE, show_violin = TRUE, 
                                 theme_size = 10) {
    ## Check object is an SCESet object
    if ( !is(object, "SCESet") )
        stop("object must be an SCESet")
    
    ## Define number of features to plot
    if (is.logical(features))
        nfeatures <- sum(features)
    else
        nfeatures <- length(features)
    
    ## Checking arguments for expression values
    use_as_exprs <- match.arg(use_as_exprs, 
                              choices = c("exprs", "norm_exprs", "counts", 
                                          "norm_counts", "tpm", "norm_tpm", 
                                          "fpkm", "norm_fpkm", "cpm", "norm_cpm"))
    if ( use_as_exprs == "counts" && is.null(counts(object)) ) {
        warning("'use_as_exprs' argument is 'counts', but counts(object) is NULL. Plotting 'exprs' values instead.")
        use_as_exprs <- "exprs"
    }
    if ( use_as_exprs == "tpm" && is.null(tpm(object)) ) {
        warning("'use_as_exprs' argument is 'tpm', but tpm(object) is NULL. Plotting 'exprs' values instead.")
        use_as_exprs <- "exprs"
    } 
    if ( use_as_exprs == "norm_tpm" && is.null(norm_tpm(object)) ) {
        warning("'use_as_exprs' argument is 'norm_tpm', but norm_tpm(object) is NULL. Plotting 'exprs' values instead.")
        use_as_exprs <- "exprs"
    } 
    if ( use_as_exprs == "fpkm" && is.null(fpkm(object)) ) {
        warning("'use_as_exprs' argument is 'fpkm', but fpkm(object) is NULL. Plotting 'exprs' values instead.")
        use_as_exprs <- "exprs"
    }
    if ( use_as_exprs == "norm_fpkm" && is.null(norm_fpkm(object)) ) {
        warning("'use_as_exprs' argument is 'norm_fpkm', but norm_fpkm(object) is NULL. Plotting 'exprs' values instead.")
        use_as_exprs <- "exprs"
    }
    if ( use_as_exprs == "cpm" && is.null(cpm(object)) ) {
        warning("'use_as_exprs' argument is 'cpm', but cpm(object) is NULL. Plotting 'exprs' values instead.")
        use_as_exprs <- "exprs"
    } 
    if ( use_as_exprs == "norm_cpm" && is.null(norm_cpm(object)) ) {
        warning("'use_as_exprs' argument is 'norm_cpm', but norm_cpm(object) is NULL. Plotting 'exprs' values instead.")
        use_as_exprs <- "exprs"
    } 
    
    ## Check arguments are valid
    if ( !(x %in% varLabels(object)) )
        stop("the argument 'x' should specify a column of pData(object) [see varLabels(object)]")
    if ( !is.null(colour_by) ) {
        if ( !(colour_by %in% varLabels(object)) )
            stop("the argument 'colour_by' should specify a column of pData(object) [see varLabels(object)]")
        if ( nlevels(as.factor(pData(object)[[colour_by]])) > 10 )
            stop("when coerced to a factor, 'colour_by' should have fewer than 10 levels")
    }
    if ( !is.null(shape_by) ) {
        if ( !(shape_by %in% varLabels(object)) )
            stop("the argument 'shape_by' should specify a column of pData(object) [see varLabels(object)]")
        if ( nlevels(as.factor(pData(object)[[shape_by]])) > 10 )
            stop("when coerced to a factor, 'shape_by' should have fewer than 10 levels")   
    }
    if ( !is.null(size_by) ) {
        if ( !(size_by %in% varLabels(object)) )
            stop("the argument 'size_by' should specify a column of pData(object) [see varLabels(object)]")
    }
    if ( typeof(features) == "character" ) {
        if ( !(all(features %in% featureNames(object))) )
            stop("when the argument 'features' is of type character, all features must be in featureNames(object)")
    }
    
    ## Define expression values to use
    if ( use_as_exprs == "exprs" || use_as_exprs == "norm_exprs" ) {
        to_melt <- switch(use_as_exprs,
                          exprs = as.matrix(exprs(object)[features, , 
                                                          drop = FALSE]),
                          norm_exprs = as.matrix(norm_exprs(object)[features, , 
                                                               drop = FALSE]))
        if ( !object@logged ) {
            to_melt <- log2(to_melt + 1)
            message("Expression values have been transformed to a log2 scale")
        }
        ylab <- switch(use_as_exprs,
                       exprs = "Expression (log2 scale)",
                       norm_exprs = "Normalised expression (log2 scale)")
    } else {
        if ( use_as_exprs == "counts" ) {
            count_mtrx <- as.matrix(counts(object))
            lib_size <- colSums(count_mtrx)
            cpm_to_plot <- log2(t(t(count_mtrx) / lib_size) + 1)
            to_melt <- as.matrix(cpm_to_plot[features, , drop = FALSE])
            ylab <- "log2(counts-per-million + 1)"
        }
        if ( use_as_exprs == "norm_counts" ) {
            count_mtrx <- as.matrix(norm_counts(object))
            lib_size <- colSums(count_mtrx)
            cpm_to_plot <- log2(t(t(count_mtrx) / lib_size) + 1)
            to_melt <- as.matrix(cpm_to_plot[features, , drop = FALSE])
            ylab <- "normalised log2(counts-per-million + 1)"
        }
        if ( use_as_exprs == "tpm" ) {
            tpm_mtrx <- as.matrix(tpm(object)[features, , drop = FALSE])
            to_melt <- log2(tpm_mtrx + 1)
            ylab <- "log2(transcripts-per-million + 1)"
        }
        if ( use_as_exprs == "norm_tpm" ) {
            tpm_mtrx <- as.matrix(norm_tpm(object)[features, , drop = FALSE])
            to_melt <- log2(tpm_mtrx + 1)
            ylab <- "normalised log2(transcripts-per-million + 1)"
        }
        if ( use_as_exprs == "fpkm" ) {
            fpkm_mtrx <- as.matrix(fpkm(object)[features, , drop = FALSE])
            to_melt <- log2(fpkm_mtrx + 1)
            ylab <- "log2(FPKM + 1)"
        }
        if ( use_as_exprs == "norm_fpkm" ) {
            fpkm_mtrx <- as.matrix(norm_fpkm(object)[features, , drop = FALSE])
            to_melt <- log2(fpkm_mtrx + 1)
            ylab <- "normalised log2(FPKM + 1)"
        }
        if ( use_as_exprs == "cpm" ) {
            cpm_mtrx <- as.matrix(cpm(object)[features, , drop = FALSE])
            to_melt <- log2(cpm_mtrx + 1)
            ylab <- "log2(counts-per-million + 1)"
        }
        if ( use_as_exprs == "norm_cpm" ) {
            cpm_mtrx <- as.matrix(norm_cpm(object)[features, , drop = FALSE])
            to_melt <- log2(cpm_mtrx + 1)
            ylab <- "normalised log2(counts-per-million + 1)"
        }
    }
    
    ## Melt the expression data and metadata into a convenient form
    evals_long <- reshape2::melt(to_melt, value.name = "evals")
    colnames(evals_long) <- c("Feature", "Cell", "evals")
    ## Extend the samples information
    samples_long <- pData(object)[rep(seq_len(ncol(object)), each = nfeatures), ]
    
    ## Construct a ggplot2 aesthetic for the plot
    aesth <- aes()
    aesth$x <- as.symbol(x)
    aesth$y <- as.symbol("evals")
    if ( !is.null(colour_by) )
        aesth$colour <- as.symbol(colour_by)
    else {
        if ( !is.null(is_exprs(object)) ) {
            ## Colour by is_exprs if we can (i.e. is_exprs is not NULL)
            isexpr_long <- reshape2::melt(is_exprs(object)[features,], 
                                          value.name = "is_exprs")
            evals_long <- dplyr::mutate(evals_long, 
                                        Is_Expressed = as.vector(isexpr_long$is_exprs))
            aesth$colour <- as.symbol("Is_Expressed")
        }
    }
    if ( !is.null(shape_by) )
        aesth$shape <- as.symbol(shape_by)
    ## Define sensible x-axis label if NULL
    if ( is.null(xlab) )
        xlab <- x
    
    ## Combine the expression values and sample information
    object <- cbind(evals_long, samples_long)
    
    ## Make the plot
    plot_out <- plotExpressionDefault(object, aesth, ncol, xlab, ylab, 
                                      show_median, show_violin)
    
    ## Define plotting theme
    if ( library(cowplot, logical.return = TRUE) )
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    else
        plot_out <- plot_out + theme_bw(theme_size)
    plot_out
}


#' @param aesth an \code{aes} object to use in the call to \code{\link{ggplot}}.
#' @param ylab character string defining a label for the y-axis (y-axes) of the 
#' plot.
#' @rdname plotExpression
#' @aliases plotExpression
#' @export
plotExpressionDefault <- function(object, aesth, ncol=2, xlab=NULL, 
                                   ylab=NULL, show_median=TRUE, show_violin=TRUE) {
    if ( !("Feature" %in% names(object)) )
        stop("object needs a column named 'Feature' to define the feature(s) by which to plot expression.")
 
    ## Define the plot
    plot_out <- ggplot(object, aesth) +
        facet_wrap(~Feature, ncol = ncol, scales = 'free_y') +
        ggthemes::scale_colour_tableau() +
        xlab(xlab) +
        ylab(ylab)
    if (show_violin) {
        plot_out <- plot_out + geom_violin(group = 1, colour = "gray80", 
                                           fill = "gray80", scale = "width")
    }
    if (show_median) {
        plot_out <- plot_out +
            stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
                         geom = "crossbar", width = 0.3, alpha = 0.8)
    }
    plot_out <- plot_out + 
        geom_jitter(size = 4, alpha = 0.8, position = position_jitter(height = 0))
    plot_out + theme_bw()
}


#' @rdname plotExpression
#' @aliases plotExpression
#' @export
setMethod("plotExpression", signature(object = "SCESet"),
          function(object, features, x, use_as_exprs = "exprs", colour_by = NULL, 
                   shape_by = NULL, size_by = NULL, ncol = 2, xlab = NULL, 
                   show_median=TRUE, show_violin=TRUE) {
              plotExpressionSCESet(object, features, x, use_as_exprs, 
                                        colour_by, shape_by, size_by, ncol, 
                                        xlab, show_median, show_violin)
          })

#' @rdname plotExpression
#' @aliases plotExpression
#' @export
setMethod("plotExpression", signature(object = "SCESet"),
          function(object, ...) {
              plotExpressionSCESet(object, ...)
          })


#' @rdname plotExpression
#' @aliases plotExpression
#' @export
setMethod("plotExpression", signature("data.frame"),
          function(object, ...) {
              plotExpressionDefault(object, ...)
          })

################################################################################

#' Plot metadata for cells or features
#'
#' @param object a data.frame (or object that can be coerced to such) object 
#' containing metadata in columns to plot.
#' @param aesth aesthetics function call to pass to ggplot. This function 
#' expects at least x and y variables to be supplied. The default is to plot 
#' coverage against log10(depth).
#' @param shape numeric scalar to define the plotting shape. Ignored if shape is
#' included in the \code{aesth} argument.
#' @param alpha numeric scalar (in the interval 0 to 1) to define the alpha 
#' level (transparency) of plotted points. Ignored if alpha is included in the 
#' \code{aesth} argument.
#' @param size numeric scalar to define the plotting size. Ignored if size is
#' included in the \code{aesth} argument.
#'
#' @details Plot cell or feature metadata from an SCESet object. If one variable
#'  is supplied then a density plot will be returned. If both variables are 
#' continuous (numeric) then a scatter plot will be returned. If one variable is
#' discrete and one continuous then a violin plot with jittered points overlaid 
#' will be returned. If both variables are discrete then a jitter plot will be 
#' produced. The object returned is a ggplot object, so further layers and 
#' plotting options (titles, facets, themes etc) can be added.
#' 
#' @export
#' 
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset <- calculateQCMetrics(example_sceset)
#' plotMetadata(pData(example_sceset))
#'
plotMetadata <- function(object, 
                         aesth=aes_string(x = "log10(depth)", y = "coverage"), 
                         shape = NULL, alpha = NULL, size = NULL) {
    ## Must have at least an x variable in the aesthetics
    if (is.null(aesth$x))
        stop("No x variable defined. Must have at least an x variable defined
             in the aesth argument.")
    
    ## Determine type of plot to make
    ### Define plot type if there is only an x variable but no y
    if (is.null(aesth$y)) {
        typeof_x <- typeof(aesth[[1]])
        plot_type <- "density"
        if (typeof_x == "symbol") {
            var_name <- as.character(aesth[[1]])
            x <- typeof(object[, var_name])
            if (is.character(x) | is.factor(x))
                plot_type <- "bar"
        }   
    } else {
        ### Define plot type if both x and y variables are provided  
        typeof_x <- typeof(aesth$x)
        typeof_y <- typeof(aesth$y)
        var_type_x <- var_type_y <- "continuous"
        if (typeof_x == "symbol") {
            var_name <- as.character(aesth$x)
            x <- object[, var_name]
            if (is.character(x) | is.factor(x))
                var_type_x <- "discrete"
        }   
        if (typeof_y == "symbol") {
            var_name <- as.character(aesth$y)
            y <- object[, var_name]
            if (is.character(y) | is.factor(y))
                var_type_x <- "discrete"
        }
        if ( var_type_x == "continuous" & var_type_y == "continuous" )
            plot_type <- "scatter"
        else {
            if ( var_type_x == "discrete" & var_type_y == "discrete" )
                plot_type <- "jitter"
            else
                plot_type <- "violin"
        }
    }
    
    ## Setup plot
    show_size_guide <- show_alpha_guide <- show_shape_guide <- TRUE
    if ( is.null(aesth$size) ) {
        show_size_guide <- FALSE
        if ( is.null(size) )
            aesth$size <- 4
        else 
            aesth$size <- size
    }
    if ( is.null(aesth$alpha) ) {
        show_alpha_guide <- FALSE
        if ( is.null(alpha) )
            aesth$alpha <- 0.7
        else
            aesth$alpha <- alpha
    }
    if ( is.null(aesth$shape) ) {
        show_shape_guide <- FALSE
        if ( is.null(shape) )
            shape <- 16
    }
    
    ## Set up basics of plot
    plot_out <- ggplot(object, aesth)
    
    ## Density plot
    if (plot_type == "bar") {
        plot_out <- plot_out + geom_bar(stat = "identity") 
    }
    if (plot_type == "density") {
        plot_out <- plot_out + geom_density(kernel = "rectangular", size = 2) +
            geom_rug(alpha = 0.5, size = 1)
    }
    if (plot_type == "jitter") {
        if ( !show_shape_guide )
            plot_out <- plot_out + geom_jitter(shape = shape)
        else 
            plot_out <- plot_out + geom_jitter()
    }
    if (plot_type == "scatter") {
        if ( !show_shape_guide )
            plot_out <- plot_out + geom_point(shape = shape)
        else
            plot_out <- plot_out + geom_point()
        plot_out <- plot_out + geom_rug(alpha = 0.5, size = 1)
    }
    if (plot_type == "violin") {
        plot_out <- plot_out + geom_violin(scale = "width")
        if ( !show_shape_guide )
            plot_out  <- plot_out + geom_jitter(shape = shape)
        else
            plot_out  <- plot_out + geom_jitter()
    }
    
    ## Define plot colours
    plot_out <- plot_out + ggthemes::scale_colour_tableau() +
        ggthemes::scale_fill_tableau()
    
    ## Define plot theme
    plot_out <- plot_out + theme_bw(16) +
        theme(legend.justification = c(0, 0), 
              legend.position = c(0, 0), 
              legend.title = element_text(size = 11),
              legend.text = element_text(size = 10))
    
    ## Tweak plot guides
#     plot_out <- plot_out + guides(colour = guide_legend(override.aes = list(size = 2)),
#            shape = guide_legend(override.aes = list(size = 2)),
#            fill = guide_legend(override.aes = list(size = 2)))
    if ( !show_alpha_guide )
    plot_out <- plot_out + guides(alpha = FALSE)
    if ( !show_shape_guide )
        plot_out <- plot_out + guides(shape = FALSE)
    if ( !show_size_guide )
        plot_out <- plot_out + guides(size = FALSE)

    ## Return plot object
    plot_out
}


################################################################################

#' Plot phenotype data from an SCESet object
#'
#' @param object an SCESet object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param aesth aesthetics function call to pass to ggplot. This function 
#' expects at least x and y variables to be supplied. The default is to plot 
#' coverage against log10(depth).
#' @param theme_size numeric scalar giving default font size for plotting theme 
#' (default is 10).
#' @param ... arguments passed to \code{\link{plotMetadata}}.
#'
#' @details Plot phenotype data from an SCESet object. If one variable is 
#' supplied then a density plot will be returned. If both variables are 
#' continuous (numeric) then a scatter plot will be returned. If one variable is
#' discrete and one continuous then a violin plot with jittered points overlaid 
#' will be returned. If both variables are discrete then a jitter plot will be 
#' produced. The object returned is a ggplot object, so further layers and 
#' plotting options (titles, facets, themes etc) can be added.
#' .
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset <- calculateQCMetrics(example_sceset)
#' plotPhenoData(example_sceset, 
#' aesth = aes_string(x = "log10(depth)", y = "coverage", colour = "Mutation_Status"))
#' 
plotPhenoData <- function(object, aesth=aes_string(x = "log10(depth)", 
                                                   y = "coverage"), 
                          theme_size = 10, ...) {
    ## We must have an SCESet object
    if (!is(object, "SCESet"))
        stop("object must be an SCESet object.")

    ## Pass pData(object) to plotMetadata
    plot_out <- plotMetadata(pData(object), aesth, ...)   
    
    ## Define plotting theme
    if ( library(cowplot, logical.return = TRUE) )
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    else
        plot_out <- plot_out + theme_bw(theme_size)
    ## Return plot object
    plot_out
}


################################################################################

#' Plot feature (gene) data from an SCESet object
#'
#' @param object an SCESet object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param aesth aesthetics function call to pass to ggplot. This function 
#' expects at least x and y variables to be supplied. The default is to produce 
#' a density plot of number of cells expressing the feature (requires 
#' \code{calculateQCMetrics} to have been run on the SCESet object prior).
#' @param theme_size numeric scalar giving default font size for plotting theme 
#' (default is 10).
#' @param ... arguments passed to \code{\link{plotMetadata}}.
#'
#' @details Plot feature (gene) data from an SCESet object. If one variable is 
#' supplied then a density plot will be returned. If both variables are 
#' continuous (numeric) then a scatter plot will be returned. If one variable is
#' discrete and one continuous then a violin plot with jittered points overlaid 
#' will be returned. If both variables are discrete then a jitter plot will be 
#' produced. The object returned is a ggplot object, so further layers and 
#' plotting options (titles, facets, themes etc) can be added.
#' 
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset <- calculateQCMetrics(example_sceset)
#' plotFeatureData(example_sceset, aesth=aes(x=n_cells_exprs, y=pct_total_counts))
#' 
plotFeatureData <- function(object, 
                            aesth = aes_string(x = "n_cells_exprs",
                                               y = "prop_total_counts"), 
                            theme_size = 10, ...) {
    ## We must have an SCESet object
    if (!is(object, "SCESet"))
        stop("object must be an SCESet object.")
    
    ## Pass pData(object) to plotMetadata
    plot_out <- plotMetadata(fData(object), aesth, ...)   

    ## Define plotting theme
    if ( library(cowplot, logical.return = TRUE) )
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    else
        plot_out <- plot_out + theme_bw(theme_size)
    ## Return plot object
    plot_out
}


################################################################################
### Multiplot function for ggplot2 plots

#' Multiple plot function for ggplot2 plots
#'
#' Place multiple \code{\link[ggplot2]{ggplot}} plots on one page. 
#'
#' @param ...,plotlist ggplot objects can be passed in ..., or to plotlist (as 
#' a list of ggplot objects)
#' @param cols numeric scalar giving the number of columns in the layout
#' @param layout a matrix specifying the layout. If present, \code{cols} is 
#' ignored.
#'
#' @details If the layout is something like 
#' \code{matrix(c(1,2,3,3), nrow=2, byrow=TRUE)}, then plot 1 will go in the 
#' upper left, 2 will go in the upper right, and 3 will go all the way across 
#' the bottom. There is no way to tweak the relative heights or widths of the 
#' plots with this simple function. It was adapted from 
#' \url{http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/}
#' 
#' @export
#' @examples 
#' library(ggplot2)
#' ## This example uses the ChickWeight dataset, which comes with ggplot2
#' ## First plot
#' p1 <- ggplot(ChickWeight, aes(x = Time, y = weight, colour = Diet, group = Chick)) +
#'    geom_line() +
#'    ggtitle("Growth curve for individual chicks")
#' ## Second plot
#' p2 <- ggplot(ChickWeight, aes(x = Time, y = weight, colour = Diet)) +
#'    geom_point(alpha = .3) +
#'    geom_smooth(alpha = .2, size = 1) +
#'    ggtitle("Fitted growth curve per diet")
#' ## Third plot
#' p3 <- ggplot(subset(ChickWeight, Time == 21), aes(x = weight, colour = Diet)) +
#'    geom_density() +
#'    ggtitle("Final weight, by diet")
#' ## Fourth plot
#' p4 <- ggplot(subset(ChickWeight, Time == 21), aes(x = weight, fill = Diet)) +
#'     geom_histogram(colour = "black", binwidth = 50) +
#'    facet_grid(Diet ~ .) +
#'    ggtitle("Final weight, by diet") +
#'    theme(legend.position = "none")        # No legend (redundant in this graph)  
#' ## Combine plots and display
#' multiplot(p1, p2, p3, p4, cols = 2)   
#' 
multiplot <- function(..., plotlist = NULL, cols = 1, layout = NULL) {
    ## Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    num_plots <- length(plots)
    
    ## If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        ## Make the panel
        ## ncol: Number of columns of plots
        ## nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(num_plots / cols)),
                         ncol = cols, nrow = ceiling(num_plots / cols))
    }
    
    if (num_plots == 1) {
        print(plots[[1]])
    } else {
        ## Set up the page
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(
            layout = grid::grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:num_plots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}










