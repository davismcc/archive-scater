## Suite of plotting functions

#' Plot expression values for a set of genes or features
#'
#' @param features a character vector of feature names or Boolean
#' vector indicating which features should have their expression
#' values plotted
#' @param data_object an SCESet object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param aesth aesthetics function call to pass to ggplot
#' @param ncol number of columns to be used for the panels of the plot
#' @param xlab label for x-axis
#' @param ylab label for y-axis
#' @param show_median logical, show the median for each group on the plot
#' @param show_violin logical, show a violin plot for the distribution
#' for each group on the plot
#'
#' @details Plot expression values (default log2(counts-per-million +
#' 1)) for a set of genes or features.
#' @export
#'
plotExpression <- function(features, data_object, aesth, ncol = 2, 
                           xlab = "Patient", 
                           ylab = "log2(counts-per-million + 1)",
                           show_median = FALSE, show_violin = FALSE) {
    ## Define number of features to plot
    if(is.logical(features))
        nfeatures <- sum(features)
    else
        nfeatures <- length(features)
    ## Melt the expression data and metadata into a convenient form
    to_melt <- as.matrix(exprs(data_object)[features, , drop = FALSE])
    evals_long <- reshape2::melt(to_melt, value.name = "evals")
    colnames(evals_long) <- c("Feature", "Cell", "evals")
    if( data_object@logged )
        log2_evals = evals_long$evals
    else
        log2_evals = log2(evals_long$evals + 1)
    evals_long <- dplyr::mutate(evals_long, log2_evals = log2_evals)
    isexpr_long <- reshape2::melt(isExpr(data_object)[features,], 
                                  value.name = "isExpr")
    evals_long <- dplyr::mutate(evals_long, 
                                IsExpressed = as.vector(isexpr_long$isExpr))
    ## Extend the samples information
    samples_long <- pData(data_object)[rep(seq_len(ncol(data_object)), 
                                            each=nfeatures), ]
    ## Combine the expression values and sample information
    object_to_plot <- cbind(evals_long, samples_long)
    ## Define the plot
    plot_out <- ggplot(object_to_plot, aesth) +
        facet_wrap(~ Feature, ncol = ncol, scales = 'free_y') +
            ggthemes::scale_colour_tableau() +
                xlab(xlab) +
                    ylab(ylab)
    if(show_violin) {
        plot_out <- plot_out + geom_violin(group = 1, colour = "gray80", 
                                           fill = "gray80")
    }
    if(show_median) {
        plot_out <- plot_out +
            stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
                         geom = "crossbar", width = 0.3, alpha = 0.8)
    }
    plot_out <- plot_out + 
        geom_jitter(size = 4, alpha = 0.8, position = position_jitter(height = 0))
    plot_out
}

