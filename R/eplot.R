## Suite of functions for plotting expression values

#' Plot expression values for a set of genes or features
#'
#' @param features a character vector of feature names or Boolean
#' vector indicating which features should have their expression
#' values plotted
#' @param data_object a DGEList object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param aes aesthetics function call to pass to ggplot
#' @param ncol number of columns to be used for the panels of the plot
#' @param show_median logical, show the median for each group on the plot
#' @param show_violin logical, show a violin plot for the distribution
#' for each group on the plot
#' 
#' @details Plot expression values (default log2(counts-per-million +
#' 1)) for a set of genes or features.
#' @export
#'
eplot <- function(features, data_object, aes, ncol = 2, xlab = "Patient", ylab = "log2(counts-per-million + 1)",
                  show_median = FALSE, show_violin = FALSE) {
    ## Define number of features to plot
    if(is.logical(features))
        nfeatures <- sum(features)
    else
        nfeatures <- length(features)
    ## Melt the expression data and metadata into a convenient form
    evals_long <- melt(data_object$cpm[features,], value.name = "evals")
    colnames(evals_long) <- c("Feature", "Cell", "evals")
    evals_long <- mutate(evals_long, log2_evals = log2(evals + 1))
    isexpr_long <- melt(data_object$isexpr[features,], value.name = "isexpr")
    evals_long <- mutate(evals_long, IsExpressed = isexpr_long$isexpr)
    ## Extend the samples information
    samples_long <- data_object$samples[rep(seq_len(nrow(data_object$samples)), each=nfeatures),]
    ## Combine the expression values and sample information
    object_to_plot <- cbind(evals_long, samples_long)
    ## Define the plot
    plot_out <- ggplot(object_to_plot, aes) +
        facet_wrap(~ Feature, ncol = ncol, scales = 'free_y') +
            scale_colour_tableau() +
                xlab(xlab) +
                    ylab(ylab)
    if(show_violin)
        plot_out <- plot_out + geom_violin(group = 1, colour = "gray80", fill = "gray80")
    if(show_median)
        plot_out <- plot_out +
            stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3, alpha = 0.8)
    plot_out <- plot_out + geom_jitter(size = 4, alpha = 0.8, position = position_jitter(height = 0))
    plot_out
}






