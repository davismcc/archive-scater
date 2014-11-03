#' Beeswarm plot for expression levels
#'
#' \code{beeplot} creates a beeswarm plot with background boxplots to
#' show the distribution of expression values for a set of genes.
#'
#' @param expr_values matrix of expression values to be plotted.
#' @param value_name label for the expression values (default:
#' "log2_cpm").
#' @param group vector or factor that defines groups to be coloured
#' disctinctly in the plot.
#' @param ylab label for the y-axis.
#' @param optional character vector indicating whether genes are up-
#' or down-regulated (or giving any other information to be plotted
#' above beeswarms).
#' @param boxplot show a boxplot behind the beeswarm? (Default = \code{TRUE}).
#' @export
#'  
bplot <- function(expr_values, value_name = "log2_cpm", group = NULL, ylab = expression("log2(cpm)"),
                  de_direction = NULL, boxplot = TRUE) {
    ## Melt dataframe for convenient ggplot2 plotting
    ntop <- nrow(expr_values)
    expr_values_long <- melt(expr_values, value.name = value_name)
    colnames(expr_values_long) <- c("Gene", "Cell", value_name)
    expr_values_long <- mutate(expr_values_long, Gene = factor(Gene, levels = rownames(expr_values)))
    ## Define beeswarm coords
    bswarm <- beeswarm::beeswarm(get(value_name) ~ Gene, data = expr_values_long, method = 'swarm',
                         corral = 'wrap', pwcol = Cell)[, c(1, 2, 4, 6)]
    colnames(bswarm) <- c("x", "y", "Cell", "Gene")
    if( !is.null(group) ) {
        if( length(group) != ncol(expr_values) )
            error(paste("group of incorrect length. Should have length", ncol(expr_values), "equal to number of cells"))
        else
            group_long <- rep(group, ntop)
        bswarm <- mutate(bswarm, Group = group_long)
    }
    ## Define DE direction labels if available
    if( !is.null(de_direction) )
    de_dir <- data.frame(xl = 1:ntop, yl = max(bswarm$y) + 1, label = de_direction)
    ## Plot beeswarm with ggplot2
    beeswarm.plot <- ggplot2::ggplot(bswarm, aes(x, y)) +
        xlab("") +
            scale_y_continuous(ylab) +
                scale_x_continuous(breaks = c(1:ntop), labels = unique(bswarm$Gene), expand = c(0, 0.5)) +
                    ggthemes::theme_few() +
                        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
    if( boxplot )
        beeswarm.plot <- beeswarm.plot +
            geom_boxplot(aes(x, y, group = round_any(x, 1, round)), outlier.shape = NA)
    if( is.null(group) )
        beeswarm.plot <- beeswarm.plot +
            geom_point(size = 3, alpha = 0.75)
    else
        beeswarm.plot <- beeswarm.plot +
            geom_point(aes(colour = Group), size = 3, alpha = 0.75)
    if( !is.null(de_direction) )
        beeswarm.plot <- beeswarm.plot +
            geom_text(aes(xl, yl), label = de_dir$label, data = de_dir)
    ## Print the plot
    print(beeswarm.plot)
    ## Return the ggplot object (which can then be saved to file)
    beeswarm.plot
}

