## Suite of plotting functions

################################################################################

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
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset <- calculateQCMetrics(example_sceset)
#' plotExpression(1:5, example_sceset, 
#' aes(x=Mutation_Status, y=log2_evals, colour=Patient), xlab="Mutation_Status",
#' show_violin=TRUE, show_median=TRUE)
#' 
plotExpression <- function(features, data_object, aesth, ncol = 2, 
                           xlab = "", 
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
    isexpr_long <- reshape2::melt(isExprs(data_object)[features,], 
                                  value.name = "isExprs")
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
    plot_out + theme_bw()
}


################################################################################

#' Plot metadata for cells or genes
#'
#' @param object a data.frame (or object that can be coerced to such) object 
#' containing metadata in columns to plot.
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
              legend.title = element_text(8),
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



