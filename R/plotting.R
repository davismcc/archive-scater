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
#'  @param y optional argument for generic \code{plot} functions, not used for 
#'  plotting an \code{SCESet} object
#'  @param ... arguments passed to \code{plot.SCESet}
#' @param ncol number of columns to use for \code{facet_wrap} if only one block is 
#' defined.
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
#' plot(example_sceset, use_as_exprs="exprs")
#' plot(example_sceset, use_as_exprs="exprs", colour_by="Cell_Cycle")
#' plot(example_sceset, use_as_exprs="exprs", block1="Treatment", 
#' colour_by="Cell_Cycle")
#' plot(example_sceset, use_as_exprs="exprs", block1="Treatment", 
#' block2="Mutation_Status", colour_by="Cell_Cycle")
#' # What happens if chosen expression values are not available?
#' plot(example_sceset, block1="Treatment", colour_by="Cell_Cycle") 
#' 
setMethod("plot", signature("SCESet"),
          function(x, ...) {
              plot.SCESet(x, ...)
          })

#' @rdname plot
#' @aliases plot
#' @export
plot.SCESet <- function(x, block1=NULL, block2=NULL, colour_by=NULL, 
                        nfeatures=500, use_as_exprs="tpm", ncol=3) {
    object <- x
    if( !is(object, "SCESet") )
        stop("Object must be an SCESet")
    if( !is.null(block1) ) {
        if( !(block1 %in% colnames(pData(object))) )
            stop("The block1 argument must either be NULL or a column of pData(object).")        
    } 
    if( !is.null(block2) ) {
        if( !(block2 %in% colnames(pData(object))) )
            stop("The block2 argument must either be NULL or a column of pData(object).")        
    }  
    if( !is.null(colour_by) ) {
        if( !(colour_by %in% colnames(pData(object))) )
            stop("The colour_by argument must either be NULL or a column of pData(object).")        
    }  
    
    ## Define an expression matrix depending on which values we're using
    use_as_exprs <- match.arg(use_as_exprs, c("exprs", "tpm", "fpkm", "counts"))
    exprs_mat <- switch(use_as_exprs,
                        exprs=exprs(object),
                        tpm=tpm(object),
                        fpkm=fpkm(object),
                        counts=counts(object))
    if( is.null(exprs_mat) ) {
        warning(paste0("The object does not contain ", use_as_exprs, " expression values. Using exprs(object) values instead."))
        exprs_mat <- exprs(object)
        use_as_exprs <- "exprs"
    }
    if( use_as_exprs=="exprs" & object@logged )
        exprs_mat <- 2^(exprs_mat)
    
    ## Use plyr to get the sequencing real estate accounted for by features
    nfeatures_total <- nrow(exprs_mat)
    seq_real_estate <- plyr::aaply(exprs_mat, 2, .fun = function(x) {
        cumsum(sort(x, decreasing=TRUE))
    }) %>% t()
    rownames(seq_real_estate) <- 1:nfeatures_total    
    nfeatures_to_plot <- nfeatures
    seq_real_estate_long <- reshape2::melt(seq_real_estate[1:nfeatures_to_plot, ], 
                                           value.name = "exprs")
    
    ## Get the proportion of the library accounted for by the top features
    prop_library <- reshape2::melt(t(t(seq_real_estate[1:nfeatures_to_plot, ]) /
                                         colSums(exprs_mat)), 
                                   value.name="prop_library")
    colnames(seq_real_estate_long) <- c("Feature", "Cell", "exprs")
    seq_real_estate_long$Proportion_Library <- prop_library$prop_library
    
    ## Add block and colour_by information if provided
    if( !is.null(block1) )
        seq_real_estate_long <- dplyr::mutate(
            seq_real_estate_long, block1=as.factor(rep(object[[block1]], 
                                                       each=nfeatures_to_plot)))
    if( !is.null(block2) )                                  
        seq_real_estate_long <- dplyr::mutate(
            seq_real_estate_long, block2=as.factor(rep(object[[block2]], 
                                                       each=nfeatures_to_plot)))
    if( !is.null(colour_by) )                                  
        seq_real_estate_long <- dplyr::mutate(
            seq_real_estate_long, colour_by=rep(object[[colour_by]], 
                                             each=nfeatures_to_plot))
    ## Set up plot
    if( is.null(colour_by) ) {
        plot_out <- ggplot(seq_real_estate_long, 
                           aes(x=Feature, y=Proportion_Library, group=Cell)) +
            geom_line(linetype="solid", alpha = 0.3)
    } else {
        plot_out <- ggplot(seq_real_estate_long, 
                           aes(x=Feature, y=Proportion_Library, group=Cell,
                               colour=colour_by)) +
            geom_line(linetype="solid", alpha = 0.3)
    }
    ## Deal with blocks for grid
    if( !(is.null(block1) | is.null(block2)) )
        plot_out <- plot_out + facet_grid(block2~block1)
    else {
        if( !is.null(block1) & is.null(block2) ) {
            plot_out <- plot_out + 
                facet_wrap(~block1, ncol=ncol)
        }
        if( is.null(block1) & !is.null(block2) ) {
            plot_out <- plot_out + 
                facet_wrap(~block2, ncol=ncol)
        }
    }
    ## Add extra plot theme and details
    plot_out <- plot_out + theme_bw(16) + 
        ggthemes::scale_colour_tableau(name=colour_by) + 
        xlab("Number of features") + ylab("Cumulative proportion of library")
    
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
#' @param ntop numeric scalar indicating the number of most variable genes to 
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
#' plotPCA(example_sceset, colour_by="Cell_Cycle")
#' plotPCA(example_sceset, colour_by="Cell_Cycle", shape_by="Treatment")
#' plotPCA(example_sceset, colour_by="Cell_Cycle", shape_by="Treatment", 
#' size_by="Mutation_Status")
#' plotPCA(example_sceset, shape_by="Treatment", size_by="Mutation_Status")
#' plotPCA(example_sceset, feature_set=1:100, colour_by="Treatment", 
#' shape_by="Mutation_Status")
#' 
#' plotPCA(example_sceset, shape_by="Treatment", return_SCESet=TRUE)
#' 
#' ## Examples plotting more than 2 PCs
#' plotPCA(example_sceset, ncomponents=8)
#' plotPCA(example_sceset, ncomponents=4, colour_by="Treatment", 
#' shape_by="Mutation_Status")
#' 
plotPCA.SCESet <- function(object, ntop=500, ncomponents=2, colour_by=NULL, 
                           shape_by=NULL, size_by=NULL, feature_set=NULL, 
                           return_SCESet=FALSE, scale_features=TRUE) {
    ## Check arguments are valid
    if( !is.null(colour_by) ) {
        if( !(colour_by %in% varLabels(object)) )
            stop("the argument 'colour_by' should specify a column of pData(object) [see varLabels(object)]")
        if( nlevels(as.factor(pData(object)[[colour_by]])) > 10 )
            stop("when coerced to a factor, 'colour_by' should have fewer than 10 levels")
    }
    if( !is.null(shape_by) ) {
        if( !(shape_by %in% varLabels(object)) )
            stop("the argument 'shape_by' should specify a column of pData(object) [see varLabels(object)]")
        if( nlevels(as.factor(pData(object)[[shape_by]])) > 10 )
            stop("when coerced to a factor, 'shape_by' should have fewer than 10 levels")   
    }
    if( !is.null(size_by) ) {
        if( !(size_by %in% varLabels(object)) )
            stop("the argument 'size_by' should specify a column of pData(object) [see varLabels(object)]")
    }
    if( !is.null(feature_set) & typeof(feature_set)=="character" ) {
        if( !(all(feature_set %in% featureNames(object))) )
            stop("when the argument 'feature_set' is of type character, all features must be in featureNames(object)")
    }
    
    ## Define features to use: either ntop, or if a set of features is defined, then those
    if( is.null(feature_set) ) {
        rv <- matrixStats::rowVars(exprs(object))
        feature_set <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
    }
    ## Standardise expression if stand_exprs(object) is null
    if( is.null(stand_exprs(object)) ) {
        stand_exprs(object) <- exprs(object) %>% t %>% 
            scale(scale=scale_features) %>% t
        message("stand_exprs(object) was null, so standardising exprs values for PCA.")
    } else
        message("stand_exprs(object) was not null, so using these values for PCA.")
    
    ## Compute PCA
    pca <- prcomp(t(stand_exprs(object)[feature_set, ]))
    percentVar <- pca$sdev^2 / sum(pca$sdev^2)
    
    ## Define data.frame for plotting
    df_to_plot <- data.frame(pca$x[, 1:ncomponents], 
                             row.names=sampleNames(object))
    if( !is.null(colour_by) )
        df_to_plot <- data.frame(df_to_plot, 
                                 colour_by=as.factor(pData(object)[[colour_by]]))
    if( !is.null(shape_by) )
        df_to_plot <- data.frame(df_to_plot, 
                                 shape_by=as.factor(pData(object)[[shape_by]]))
    if( !is.null(size_by) )
        df_to_plot <- data.frame(df_to_plot, size_by=pData(object)[[size_by]])
    
    plot_out <- plotReducedDim.default(df_to_plot, ncomponents, colour_by, shape_by, 
                                       size_by, percentVar)
    
    ## Plot PCA and return appropriate object
    if(return_SCESet) {
        df_out <- pca$x[, 1:ncomponents]
        rownames(df_out) <- sampleNames(object)
        attr(df_out, "percentVar") <- percentVar[1:ncomponents]
        reducedDimension(object) <- df_out
        print(plot_out)
        return(object)
    } else {
        ## Return PCA plot
        plot_out   
    }
}
    

#' @rdname plotPCA
#' @aliases plotPCA
#' @export
setMethod("plotPCA", signature("SCESet"),
          function(object, ntop=500, ncomponents=2, colour_by=NULL, shape_by=NULL, 
                   size_by=NULL, feature_set=NULL, return_SCESet=FALSE, 
                   scale_features=TRUE) {
              plotPCA.SCESet(object, ntop, ncomponents, colour_by, shape_by, size_by, 
                             feature_set, return_SCESet, scale_features)
          })


.makePairs <- function(data_matrix) {
    ## with thanks to Gaston Sanchez, who posted this code online
    ## https://gastonsanchez.wordpress.com/2012/08/27/scatterplot-matrices-with-ggplot/
    grid <- expand.grid(x=1:ncol(data_matrix), y=1:ncol(data_matrix))
    grid <- subset(grid, x != y)
    all_panels <- do.call("rbind", lapply(1:nrow(grid), function(i) {
        xcol <- grid[i, "x"]
        ycol <- grid[i, "y"]
        data.frame(xvar=names(data_matrix)[ycol], yvar=names(data_matrix)[xcol], 
                   x=data_matrix[, xcol], y=data_matrix[, ycol], data_matrix)
    }))
    all_panels$xvar <- factor(all_panels$xvar, levels=names(data_matrix))
    all_panels$yvar <- factor(all_panels$yvar, levels=names(data_matrix))
    densities <- do.call("rbind", lapply(1:ncol(data_matrix), function(i) {
        data.frame(xvar=names(data_matrix)[i], yvar=names(data_matrix)[i], 
                   x=data_matrix[, i])
    }))
    list(all=all_panels, densities=densities)
}


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
                           shape_by=NULL, size_by=NULL, percentVar=NULL) {
    ## Define plot
    if( ncomponents > 2 ) {
        ## expanding numeric columns for pairs plot
        df_to_expand <- df_to_plot[, 1:ncomponents]
        if( is.null(percentVar) ) {
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
        plot_out <- ggplot(df_to_plot_big, aes_string(x="x", y="y")) + 
            facet_grid(xvar ~ yvar, scales="free") + 
            stat_density(aes(x=x, y=(..scaled.. * diff(range(x)) + min(x))), 
                         data=gg1$densities, position="identity", 
                         colour="grey20", geom="line") +
            theme_bw(14)
    } else {
        comps <- colnames(df_to_plot)[1:2]
        if( is.null(percentVar) ) {
            x_lab <- "Dimension 1"
            y_lab <- "Dimension 2"
        } else {
            x_lab <- paste0("Component 1: ", round(percentVar[1] * 100), "% variance")
            y_lab <- paste0("Component 2: ", round(percentVar[2] * 100), "% variance")
        }
        plot_out <- ggplot(data=df_to_plot, aes_string(x=comps[1], y=comps[2])) + 
            xlab(x_lab) + 
            ylab(y_lab) +
            geom_rug(colour="gray20", alpha=0.5) +
            theme_bw(14)        
    }
    
    ## Apply colour_by, shape_by and size_by variables if defined
    if( !is.null(colour_by) & !is.null(shape_by) & !is.null(size_by) ) {
        plot_out <- plot_out + 
            geom_point(aes_string(colour="colour_by", shape="shape_by", size="size_by"),
                       alpha=0.5) +
            ggthemes::scale_colour_tableau(name=colour_by) +
            guides(size=guide_legend(title=size_by), shape=guide_legend(title=shape_by))
    } else {
        if( sum(is.null(colour_by) + is.null(shape_by) + is.null(size_by)) == 1 ) {
            if( !is.null(colour_by) & !is.null(shape_by) ) {
                plot_out <- plot_out + 
                    geom_point(aes_string(colour="colour_by", shape="shape_by"), size=4,
                               alpha=0.5) +
                    ggthemes::scale_colour_tableau(name=colour_by) +
                    guides(shape=guide_legend(title=shape_by))
            }
            if( !is.null(colour_by) & !is.null(size_by) ) {
                plot_out <- plot_out + 
                    geom_point(aes_string(fill="colour_by", size="size_by"), 
                               shape=21, colour="gray70", alpha=0.5) +
                    ggthemes::scale_fill_tableau(name=colour_by) +
                    guides(size=guide_legend(title=size_by))
            }
            if( !is.null(shape_by) & !is.null(size_by) ) {
                plot_out <- plot_out + 
                    geom_point(aes_string(shape="shape_by", size="size_by"), 
                               fill="gray20", colour="gray20", alpha=0.5) +
                    guides(size=guide_legend(title=size_by), shape=guide_legend(title=shape_by))
            } 
        } else {
            if( sum(is.null(colour_by) + is.null(shape_by) + is.null(size_by)) == 2 ) {
                if( !is.null(colour_by) ) {
                    plot_out <- plot_out + 
                        geom_point(aes_string(fill="colour_by"), size=4, shape=21, 
                                   colour="gray70", alpha=0.5) +
                        ggthemes::scale_fill_tableau(name=colour_by)
                }
                if( !is.null(shape_by) ) {
                    plot_out <- plot_out + 
                        geom_point(aes_string(shape="shape_by"), size=4, 
                                   colour="gray20", alpha=0.5) +
                        guides(shape=guide_legend(title=shape_by))
                }
                if( !is.null(size_by) ) {
                    plot_out <- plot_out + 
                        geom_point(aes_string(size="size_by"), fill="gray20", shape=21, 
                                   colour="gray70", alpha=0.5) +
                        guides(size=guide_legend(title=size_by))
                } 
            } else {
                plot_out <- plot_out + 
                    geom_point(size=4, fill="gray20", shape=21, colour="gray70",
                               alpha=0.5)
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
                                  shape_by=NULL, size_by=NULL) {
    ## Check arguments are valid
    if( !is.null(colour_by) ) {
        if( !(colour_by %in% varLabels(object)) )
            stop("the argument 'colour_by' should specify a column of pData(object) [see varLabels(object)]")
        if( nlevels(as.factor(pData(object)[[colour_by]])) > 10 )
            stop("when coerced to a factor, 'colour_by' should have fewer than 10 levels")
    }
    if( !is.null(shape_by) ) {
        if( !(shape_by %in% varLabels(object)) )
            stop("the argument 'shape_by' should specify a column of pData(object) [see varLabels(object)]")
        if( nlevels(as.factor(pData(object)[[shape_by]])) > 10 )
            stop("when coerced to a factor, 'shape_by' should have fewer than 10 levels")   
    }
    if( !is.null(size_by) ) {
        if( !(size_by %in% varLabels(object)) )
            stop("the argument 'size_by' should specify a column of pData(object) [see varLabels(object)]")
    }
    
    ## Extract reduced dimension representation of cells
    if( is.null(reducedDimension(object)) )
        stop("reducedDimension slot of object is NULL. Need non null reducedDimension to plot.")
    red_dim <- redDim(object)
    if( ncomponents > ncol(red_dim) )
        stop("ncomponents to plot is larger than number of columns of reducedDimension(object)")
    
    ## Define data.frame for plotting
    df_to_plot <- data.frame(red_dim[, 1:ncomponents])
    if( !is.null(colour_by) )
        df_to_plot <- data.frame(df_to_plot, 
                                 colour_by=as.factor(pData(object)[[colour_by]]))
    if( !is.null(shape_by) )
        df_to_plot <- data.frame(df_to_plot, 
                                 shape_by=as.factor(pData(object)[[shape_by]]))
    if( !is.null(size_by) )
        df_to_plot <- data.frame(df_to_plot, size_by=pData(object)[[size_by]])
    
    ## Call default method to make the plot
    plotReducedDim.default(df_to_plot, ncomponents, colour_by, shape_by, size_by,
                                 percentVar=NULL)
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
#' \code{plotExpression.default} method, the \code{object} argument is a 
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
#' @param ... optional arguments (from those listed above) passed to 
#' \code{plotExpression.SCESet} or \code{plotExpression.default}
#'
#' @details Plot expression values (default log2(counts-per-million +
#' 1)) for a set of genes or features.
#' 
#' @name plotExpression
#' @aliases plotExpression plotExpression,SCESet-method plotExpression,data.frame-method
#' @import ggplot2
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
plotExpression.SCESet <- function(object, features, x, use_as_exprs="exprs", 
                                  colour_by=NULL, shape_by=NULL, size_by=NULL, 
                                  ncol=2, xlab=NULL, show_median=FALSE, 
                                  show_violin=FALSE) {
    ## Check object is an SCESet object
    if( !is(object, "SCESet") )
        stop("object must be an SCESet")
    
    ## Define number of features to plot
    if(is.logical(features))
        nfeatures <- sum(features)
    else
        nfeatures <- length(features)
    
    ## Checking arguments for expresion values
    use_as_exprs <- match.arg(use_as_exprs, 
                              choices=c("exprs", "counts", "tpm", "fpkm"))
    if( !(use_as_exprs == "exprs" | use_as_exprs == "counts" | use_as_exprs == "fpkm" | use_as_exprs == "fpkm") )
        stop("The argument 'use_as_exprs' must either be 'tpm' (transcripts-per-million), 'fpkm', exprs' to use expression values, or 'counts', to use read counts.")
    if( use_as_exprs == "counts" & is.null(counts(object)) ) {
        warning("'use_as_exprs' argument is 'counts', but counts(object) is NULL. Plotting 'exprs' values instead.")
        use_as_exprs <- "exprs"
    }
    if( use_as_exprs == "tpm" & is.null(tpm(object)) ) {
        warning("'use_as_exprs' argument is 'tpm', but tpm(object) is NULL. Plotting 'exprs' values instead.")
        use_as_exprs <- "exprs"
    } 
    if( use_as_exprs == "fpkm" & is.null(fpkm(object)) ) {
        warning("'use_as_exprs' argument is 'fpkm', but fpkm(object) is NULL. Plotting 'exprs' values instead.")
        use_as_exprs <- "exprs"
    }
    
    ## Check arguments are valid
    if( !(x %in% varLabels(object)) )
        stop("the argument 'x' should specify a column of pData(object) [see varLabels(object)]")
    if( !is.null(colour_by) ) {
        if( !(colour_by %in% varLabels(object)) )
            stop("the argument 'colour_by' should specify a column of pData(object) [see varLabels(object)]")
        if( nlevels(as.factor(pData(object)[[colour_by]])) > 10 )
            stop("when coerced to a factor, 'colour_by' should have fewer than 10 levels")
    }
    if( !is.null(shape_by) ) {
        if( !(shape_by %in% varLabels(object)) )
            stop("the argument 'shape_by' should specify a column of pData(object) [see varLabels(object)]")
        if( nlevels(as.factor(pData(object)[[shape_by]])) > 10 )
            stop("when coerced to a factor, 'shape_by' should have fewer than 10 levels")   
    }
    if( !is.null(size_by) ) {
        if( !(size_by %in% varLabels(object)) )
            stop("the argument 'size_by' should specify a column of pData(object) [see varLabels(object)]")
    }
    if( typeof(features) == "character" ) {
        if( !(all(features %in% featureNames(object))) )
            stop("when the argument 'features' is of type character, all features must be in featureNames(object)")
    }
    
    ## Define expression values to use
    if( use_as_exprs == "exprs" ) {
        to_melt <- as.matrix(exprs(object)[features, , drop = FALSE])
        if( !object@logged ) {
            to_melt <- log2(to_melt + 1)
            message("Expression values have been transformed to a log2 scale")
        }
        ylab <- "Expression (log2 scale)"
    }
    else {
        if( use_as_exprs == "counts" ) {
            count_mtrx <- as.matrix(counts(object))
            lib_size <- colSums(count_mtrx)
            cpm_to_plot <- log2(t(t(count_mtrx) / lib_size) + 1)
            to_melt <- as.matrix(cpm_to_plot[features, , drop = FALSE])
            ylab <- "log2(counts-per-million + 1)"
        }
        if( use_as_exprs == "tpm" ) {
            tpm_mtrx <- as.matrix(tpm(object)[features, , drop = FALSE])
            to_melt <- log2(tpm_mtrx + 1)
            ylab <- "log2(transcripts-per-million + 1)"
        }
        if( use_as_exprs == "fpkm" ) {
            fpkm_mtrx <- as.matrix(fpkm(object)[features, , drop = FALSE])
            to_melt <- log2(fpkm_mtrx + 1)
            ylab <- "log2(FPKM + 1)"
        }
    }
    
    ## Melt the expression data and metadata into a convenient form
    evals_long <- reshape2::melt(to_melt, value.name="evals")
    colnames(evals_long) <- c("Feature", "Cell", "evals")
    ## Extend the samples information
    samples_long <- pData(object)[rep(seq_len(ncol(object)), each=nfeatures), ]
    
    ## Construct a ggplot2 aesthetic for the plot
    aesth <- aes()
    aesth$x <- as.symbol(x)
    aesth$y <- as.symbol("evals")
    if( !is.null(colour_by) )
        aesth$colour <- as.symbol(colour_by)
    else {
        if( !is.null(is_exprs(object)) ) {
            ## Colour by is_exprs if we can (i.e. is_exprs is not NULL)
            isexpr_long <- reshape2::melt(is_exprs(object)[features,], 
                                          value.name = "is_exprs")
            evals_long <- dplyr::mutate(evals_long, 
                                        Is_Expressed=as.vector(isexpr_long$is_exprs))
            aesth$colour <- as.symbol("Is_Expressed")
        }
    }
    if( !is.null(shape_by) )
        aesth$shape <- as.symbol(shape_by)
    ## Define sensible x-axis label if NULL
    if( is.null(xlab) )
        xlab <- x
    
    ## Combine the expression values and sample information
    object <- cbind(evals_long, samples_long)
    
    ## Make the plot
    plotExpression.default(object, aesth, ncol, xlab, ylab, show_median, show_violin)
}


#' @param aesth an \code{aes} object to use in the call to \code{\link{ggplot}}.
#' @param ylab character string defining a label for the y-axis (y-axes) of the 
#' plot.
#' @rdname plotExpression
#' @aliases plotExpression
#' @export
plotExpression.default <- function(object, aesth, ncol=2, xlab=NULL, 
                                   ylab=NULL, show_median=FALSE, show_violin=FALSE) {
    if( !("Feature" %in% names(object)) )
        stop("object needs a column named 'Feature' to define the feature(s) by which to plot expression.")
 
    ## Define the plot
    plot_out <- ggplot(object, aesth) +
        facet_wrap(~Feature, ncol=ncol, scales='free_y') +
        ggthemes::scale_colour_tableau() +
        xlab(xlab) +
        ylab(ylab)
    if(show_violin) {
        plot_out <- plot_out + geom_violin(group=1, colour="gray80", 
                                           fill="gray80")
    }
    if(show_median) {
        plot_out <- plot_out +
            stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, 
                         geom="crossbar", width=0.3, alpha=0.8)
    }
    plot_out <- plot_out + 
        geom_jitter(size=4, alpha=0.8, position=position_jitter(height=0))
    plot_out + theme_bw()
}


#' @rdname plotExpression
#' @aliases plotExpression
#' @export
setMethod("plotExpression", signature(object="SCESet"),
          function(object, features, x, use_as_exprs="exprs", colour_by=NULL, 
                   shape_by=NULL, size_by=NULL, ncol=2, xlab=NULL, 
                   show_median=FALSE, show_violin=FALSE) {
              plotExpression.SCESet(object, features, x, use_as_exprs, 
                                        colour_by, shape_by, size_by, ncol, 
                                        xlab, show_median, show_violin)
          })

#' @rdname plotExpression
#' @aliases plotExpression
#' @export
setMethod("plotExpression", signature(object="SCESet"),
          function(object, ...) {
              plotExpression.SCESet(object, ...)
          })


#' @rdname plotExpression
#' @aliases plotExpression
#' @export
setMethod("plotExpression", signature("data.frame"),
          function(object, ...) {
              plotExpression.default(object, ...)
          })

################################################################################

#' Plot metadata for cells or genes
#'
#' @param object a data.frame (or object that can be coerced to such) object 
#' containing metadata in columns to plot.
#' @param aesth aesthetics function call to pass to ggplot. This function 
#' expects at least x and y variables to be supplied. The default is to plot 
#' coverage against log10(depth).
#'
#' @details Plot cell or feature metadata from an SCESet object. If one variable
#'  is supplied then a density plot will be returned. If both variables are 
#' continuous (numeric) then a scatter plot will be returned. If one variable is
#' discrete and one continuous then a violin plot with jittered points overlaid 
#' will be returned. If both variables are discrete then a jitter plot will be 
#' produced. The object returned is a ggplot object, so further layers and 
#' plotting options (titles, facets, themes etc) can be added.
#' 
#' @import magrittr
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
plotMetadata <- function(object, aesth=aes(x=log10(depth), y=coverage)) {
    ## Must have at least an x variable in the aesthetics
    if(is.null(aesth$x))
        stop("No x variable defined. Must have at least an x variable defined
             in the aesth argument.")
    
    ## Determine type of plot to make
    ### Define plot type if there is only an x variable but no y
    if(is.null(aesth$y)) {
        typeof_x <- aesth[[1]] %>% typeof
        plot_type <- "density"
        if(typeof_x == "symbol") {
            var_name <- aesth[[1]] %>% as.character
            x <- object[, var_name] %>% typeof
            if(is.character(x) | is.factor(x))
                plot_type <- "bar"
        }   
    } else {
        ### Define plot type if both x and y variables are provided  
        typeof_x <- aesth$x %>% typeof
        typeof_y <- aesth$y %>% typeof
        var_type_x <- var_type_y <- "continuous"
        if(typeof_x == "symbol") {
            var_name <- aesth$x %>% as.character
            x <- object[, var_name]
            if(is.character(x) | is.factor(x))
                var_type_x <- "discrete"
        }   
        if(typeof_y == "symbol") {
            var_name <- aesth$y %>% as.character
            y <- object[, var_name]
            if(is.character(y) | is.factor(y))
                var_type_x <- "discrete"
        }
        if( var_type_x == "continuous" & var_type_y == "continuous" )
            plot_type <- "scatter"
        else {
            if( var_type_x == "discrete" & var_type_y == "discrete" )
                plot_type <- "jitter"
            else
                plot_type <- "violin"
        }
    }
    
    ## Setup plot
    plot_out <- ggplot(object, aesth)
       
    ## Density plot
    if(plot_type == "bar") {
        plot_out <- plot_out + geom_bar(stat="identity") 
    }
    if(plot_type == "density") {
        plot_out <- plot_out + geom_density(kernel = "rectangular", size=2) +
            geom_rug(alpha=0.5)
    }
    if(plot_type == "jitter") {
        plot_out <- plot_out + geom_jitter(size=4, alpha=0.7)
    }
    if(plot_type == "scatter") {
        plot_out <- plot_out + 
            geom_point(size=5, alpha=0.7) +
            geom_rug(alpha=0.5)
            
    }
    if(plot_type == "violin") {
        plot_out <- plot_out + 
            geom_violin() +
            geom_jitter(size=4, alpha=0.7)
    }
    
    ## Define plot colours
    plot_out <- plot_out + ggthemes::scale_colour_tableau() +
        ggthemes::scale_fill_tableau()
    
    ## Define plot theme
    plot_out <- plot_out + theme_bw(16) +
        theme(legend.justification=c(0, 0), 
              legend.position=c(0, 0), 
              legend.title = element_text(size=8),
              legend.text=element_text(size=7))
    
    ## Tweak plot guides
    plot_out <- plot_out + guides(colour=guide_legend(override.aes=list(size=2)),
           shape=guide_legend(override.aes=list(size=2)),
           fill=guide_legend(override.aes=list(size=2)))
    
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
#' #plotPhenoData(example_sceset, aesth=aes(x=log10(depth), y=coverage, colour=Mutation_Status))
#' 
plotPhenoData <- function(object, aesth=aes(x=log10(depth), y=coverage)) {
    ## We must have an SCESet object
    if(!is(object, "SCESet"))
        stop("object must be an SCESet object.")

    ## Pass pData(object) to plotMetadata
    plot_out <- plotMetadata(pData(object), aesth)   
    
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
#' a density plot of number of cells expressing the gene (requires 
#' \code{calculateQCMetrics} to have been run on the SCESet object prior).
#'
#' @details Plot feature (gene) data from an SCESet object. If one variable is 
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
#' plotFeatureData(example_sceset, aesth=aes(x=n_cells_exprs, y=prop_total_reads))
#' 
plotFeatureData <- function(object, aesth=aes(x=n_cells_exprs, y=prop_total_reads)) {
    ## We must have an SCESet object
    if(!is(object, "SCESet"))
        stop("object must be an SCESet object.")
    
    ## Pass pData(object) to plotMetadata
    plot_out <- plotMetadata(fData(object), aesth)   
    
    ## Return plot object
    plot_out
}



