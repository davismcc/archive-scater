### Not currently used


#' Beeswarm plot for expression levels
#'
#' \code{bplot} creates a beeswarm plot with background boxplots to
#' show the distribution of expression values for a set of genes.
#'
#' @param expr_values matrix of expression values to be plotted.
#' @param group vector or factor that defines groups to be coloured
#' disctinctly in the plot.
#' @param ylab label for the y-axis.
#' @param de_direction optional character vector indicating whether genes are up-
#' or down-regulated (or giving any other information to be plotted
#' above beeswarms).
#' @param distribution_plot show either a violinplot or a boxplot
#' behind the beeswarm? (Default = \code{"violinplot"}).
#' @export
#'
# plotBeeswarm <- function(expr_values, group = NULL, ylab = expression("log2(cpm)"),
#                   de_direction = NULL, distribution_plot = "violinplot", col_theme = "few") {
#     library(beeswarm)
#     library(ggplot2)
#     library(ggthemes)
#     ## Melt dataframe for convenient ggplot2 plotting
#     ntop <- nrow(expr_values)
#     expr_values_long <- melt(expr_values, value.name = "Expression")
#     colnames(expr_values_long) <- c("Gene", "Cell", "Expression")
#     ## Define group factor to colour points
#     if( !is.null(group) ) {
#         if( length(group) == ncol(expr_values) )
#             group_long <- rep(group, ntop)
#         else if( length(group) == length(expr_values) )
#             group_long <- as.vector(group)
#         else
#             stop(paste("group of incorrect length. Should have length", ncol(expr_values), "equal to number of cells"))
#     } else
#         group_long <- rep(1, length(expr_values))
#     ## Add gene and group factors to expression data frame
#     expr_values_long <- mutate(expr_values_long, Gene = factor(Gene, levels = rownames(expr_values)),
#                                Group = factor(group_long))
#     expr_values_long <- arrange(expr_values_long, Gene)
#     ## Define beeswarm coords
#     bswarm <- beeswarm(Expression ~ Gene, data = expr_values_long, method = 'swarm', corral = 'wrap',
#                        corralWidth = 0.8, pwcol = Group, do.plot = FALSE)[, c(1, 2, 4, 6)]
#     colnames(bswarm) <- c("x", "y", "Group", "Gene")
#     ## Define DE direction labels if available
#     if( !is.null(de_direction) )
#         de_dir <- data.frame(xl = 1:ntop, yl = max(bswarm$y) + 1, label = de_direction)
#     ## Plot beeswarm with ggplot2
#     ## beeswarm.plot <- ggplot2::ggplot(expr_values_long, aes(x = Gene, y = Expression)) +
#     beeswarm.plot <- ggplot(bswarm, aes(x, y, colour = Group, fill = Group)) +
#         xlab("") +
#         scale_y_continuous(ylab) +
#         ## scale_x_continuous(breaks = c(1:ntop), labels = unique(bswarm$Gene), expand = c(0, 0.5)) +
#         theme_few() +
#         theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
#     if( distribution_plot == "boxplot" )
#         beeswarm.plot <- beeswarm.plot +
#         geom_boxplot(aes(x, y, group = round_any(x, 1, round)), outlier.shape = NA)
#     if( distribution_plot == "violinplot" )
#         beeswarm.plot <- beeswarm.plot +
#         geom_violin(aes(x = Gene, y = Expression, colour = Group, fill = Group),
#                     outlier.shape = NA, data = expr_values_long, trim = TRUE, position = "identity") +
#         stat_summary(aes(x = Gene, y = Expression, colour = Group),
#                      fun.y = median, fun.ymin = median, fun.ymax = median,
#                      geom="crossbar", width = 0.6, data = expr_values_long) +
#         scale_fill_manual(values = gray.colors(nlevels(expr_values_long$Group),
#                                                start = 0.6), guide = FALSE)
#     if( is.null(group) )
#         beeswarm.plot <- beeswarm.plot +
#         geom_point(size = 3, alpha = 0.75)
#     else
#         beeswarm.plot <- beeswarm.plot +
#         geom_point(size = 3, alpha = 0.75)
#     if( !is.null(de_direction) )
#         beeswarm.plot <- beeswarm.plot +
#         geom_text(aes(xl, yl), label = de_dir$label, data = de_dir)
#     ## Select the colour theme
#     switch(col_theme,
#            few = {beeswarm.plot <- beeswarm.plot + scale_colour_few()},
#            solarized_light = {beeswarm.plot <- beeswarm.plot + scale_colour_solarized("blue")},
#            solarized_dark = {beeswarm.plot <- beeswarm.plot + scale_colour_solarized("red")},
#            pander = {beeswarm.plot <- beeswarm.plot + scale_colour_pander()},
#            tableau = {beeswarm.plot <- beeswarm.plot + scale_colour_tableau("colorblind10")},
#            economist = {beeswarm.plot <- beeswarm.plot + scale_colour_economist()},
#            wsj = {beeswarm.plot <- beeswarm.plot + scale_colour_wsj("colors6", "")},
# {beeswarm.plot <- beeswarm.plot + scale_colour_few()}
#     )
# 
# ## Print the plot
# print(beeswarm.plot)
# ## Return the ggplot object (which can then be saved to file)
# beeswarm.plot
# }
# 
