## Convenience function for computing QC metrics and adding to pData & fData
### This file contains definitions for the following functions:
### * calculateQCMetrics
### * findImportantPCs
### * plotExplanatoryVariables
### * plotHighestReadCounts
### * plotQC
### 
### * .calculateSilhouetteWidth
### * .getRSquared
### * .getTypeOfVariable

#### Some ideas for adding to QC metrics - proportion of library from top x genes
#     subset(seq_real_estate_long, Gene == 10) %>% 
#         ggplot(aes(x = Proportion_Library, colour = culture)) + 
#         geom_density(aes(fill = culture), alpha = 0.5) + 
#         facet_wrap(~perturbed, ncol=2) + 
#         theme_igray(16) + scale_colour_tableau() + scale_fill_tableau() + 
#         xlab("Cumulative proportion of library from 10 most expressed genes") + 
#         ylab("Density")

#     subset(seq_real_estate_long, Gene == 50) %>% 
#         ggplot(aes(x = norm.lib.size, y = Proportion_Library, 
#                    colour = culture)) + geom_point(size = 4, alpha = 0.7) + geom_rug(alpha = 0.7) + 
#         facet_wrap(~perturbed, ncol = 2) + theme_igray(16) + scale_colour_tableau() + 
#         scale_fill_tableau() + ylab("Cumulative proportion of library from 50 most expressed genes") + 
#         xlab("Normalised library size")

#' Calculate QC metrics
#'
#' @param object an SCESet object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param gene_controls a character vector of gene names, or a logical vector, 
#' or a numeric vector of indices used to identify gene controls (for example,
#' ERCC spike-in genes, mitochondrial genes, etc).
#' @param cell_controls a character vector of cell (sample) names, or a logical 
#' vector, or a numeric vector of indices used to identify cell controls (for 
#' example, blank wells or bulk controls).
#' @param nmads numeric scalar giving the number of median absolute deviations to be 
#' used to flag potentially problematic cells based on depth (total number of 
#' counts for the cell, or library size) and coverage (number of genes with 
#' non-zero expression). Default value is 5.
#'
#' @details Calculate useful quality control metrics to help with pre-processing
#' of data and identification of potentially problematic genes and cells.
#' @importFrom Biobase pData
#' @importFrom Biobase fData
#' @importFrom Biobase exprs
#' @importFrom Biobase sampleNames<-
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data=sc_example_cell_info)
#' rownames(pd) <- pd$Cell
#' example_sceset <- newSCESet(countData=sc_example_counts, phenoData=pd)
#' example_sceset <- calculateQCMetrics(example_sceset)
#' 
calculateQCMetrics <- function(object, gene_controls=NULL, cell_controls=NULL, 
                               nmads=5) {
    ## We must have an SCESet object
    if(!is(object, "SCESet"))
        stop("object must be an SCESet object.")
    
    ## Compute cell-level metrics
    if( is.null(is_exprs(object)) ) {
        if(is.null(counts(object))) {
            stop("Please define is_exprs(object). E.g. use is_exprs(object) <- exprs(object) > 0.1")
        } else {            
            warning("is_exprs(object) is null. Defining 'is_exprs' using count data and 
a lower count threshold of 0.")   
            isexprs <- calcIsExprs(object, lowerDetectionLimit=0)
            rownames(isexprs) <- rownames(counts(object))
            colnames(isexprs) <- colnames(counts(object))
            is_exprs(object) <- isexprs
        }       
    }    
    ## Compute depth and coverage and find outliers
    depth <- colSums(counts(object))
    coverage <- colSums(is_exprs(object))
    mad_coverage <- mad(coverage)
    med_coverage <- median(coverage)
    keep_coverage <- c(med_coverage - 5*mad_coverage, 
                       med_coverage + 5*mad_coverage)
    mad_depth <- mad(log10(depth))
    med_depth <- median(log10(depth))
    keep_depth <- c(med_depth - 5*mad_depth, med_depth + 5*mad_depth)
    filter_on_depth <- (findInterval(log10(depth), keep_depth) != 1)
    filter_on_coverage <- (findInterval(coverage, keep_coverage) != 1)
    
    ## Contributions from control genes
    ### Determine if vector or list
    if( is.null(gene_controls) | length(gene_controls) == 0 ) {
        reads_from_gene_controls <- pct_reads_from_gene_controls <- 
            rep(0, ncol(object))
        is_gene_control <- rep(FALSE, nrow(object))    
        gene_controls_pdata <- data.frame(reads_from_gene_controls,
                                          pct_reads_from_gene_controls) 
        gene_controls_fdata <- data.frame(is_gene_control)
        n_sets_gene_controls <- 1
    } else {
        if( is.list(gene_controls) ) {
            gene_controls_list <- gene_controls
            n_sets_gene_controls <- length(gene_controls)
        }
        else {
            gene_controls_list <- list(gene_controls)
            n_sets_gene_controls <- 1
        }
        ## Cycle through the gene_controls list and add QC info
        for( i in seq_len(length(gene_controls_list)) ) {
            gc_set <- gene_controls_list[[i]]
            set_name <- names(gene_controls_list)[i]
            if( is.logical(gc_set) ) {
                is_gene_control <- gc_set
                gc_set <- which(gc_set)
            } else {
                is_gene_control <- rep(FALSE, nrow(object))    
            }
            if(is.character(gc_set))
                gc_set <- which(rownames(object) %in% gc_set)
            reads_from_gene_controls <- colSums(counts(object)[gc_set,])
            is_gene_control[gc_set] <- TRUE
            pct_reads_from_gene_controls <- 100 * reads_from_gene_controls / depth
            ## Construct data.frame for pData from this gene control set
            df_pdata_this <- data.frame(reads_from_gene_controls,
                                        pct_reads_from_gene_controls)
            colnames(df_pdata_this) <- paste(colnames(df_pdata_this), set_name, 
                                             sep="_")
            if( exists("gene_controls_pdata") )
                gene_controls_pdata <- cbind(gene_controls_pdata, df_pdata_this)
            else
                gene_controls_pdata <- df_pdata_this
            ## Construct data.frame for fData from this gene control set
            df_fdata_this <- data.frame(is_gene_control)
            colnames(df_fdata_this) <- paste(colnames(df_fdata_this), set_name, 
                                             sep="_")
            if( exists("gene_controls_fdata") )
                gene_controls_fdata <- cbind(gene_controls_fdata, df_fdata_this)
            else 
                gene_controls_fdata <- df_fdata_this
        }
    }
    
    ## Fix column names and define gene controls across all control sets
    if( n_sets_gene_controls == 1 ) {
        colnames(gene_controls_pdata) <- c("reads_from_gene_controls",
                                           "pct_reads_from_gene_controls")
        colnames(gene_controls_fdata) <- "is_gene_control"
    } else {
        ## Combine all gene controls
        is_gene_control <- apply(gene_controls_fdata, 1, any)
        gene_controls_fdata <- cbind(gene_controls_fdata, is_gene_control)
        ## Compute metrics using all gene controls
        reads_from_gene_controls <- colSums(counts(object)[is_gene_control,])
        pct_reads_from_gene_controls <- 100 * reads_from_gene_controls / depth
        gene_controls_pdata <- cbind(gene_controls_pdata, 
                                     reads_from_gene_controls,
                                     pct_reads_from_gene_controls)
    }
    
    ## Define log10 read counts from gene controls
    read_counts_columns <- grep("^reads_from_gene", 
                                colnames(gene_controls_pdata))
    log10_reads_from_gene_controls <- log10(
        gene_controls_pdata[, read_counts_columns, drop=FALSE] + 1)
    colnames(log10_reads_from_gene_controls) <- 
        paste0("log10_", colnames(gene_controls_pdata)[read_counts_columns])
    ## Define reads from biological genes
    reads_from_biological_genes <- depth - 
        gene_controls_pdata$reads_from_gene_controls
    
    ## Define cell controls
    ### Determine if vector or list
    if( is.null(cell_controls) | length(cell_controls) == 0 ) {
        is_cell_control <- rep(FALSE, ncol(object))     
        cell_controls_pdata <- data.frame(is_cell_control)
        n_sets_cell_controls <- 1
    } else {
        if( is.list(cell_controls) ) {
            cell_controls_list <- cell_controls
            n_sets_cell_controls <- length(cell_controls)
        }
        else {
            cell_controls_list <- list(cell_controls)
            n_sets_cell_controls <- 1
        }
        for( i in seq_len(n_sets_cell_controls) ) {
            cc_set <- cell_controls_list[[i]]
            set_name <- names(cell_controls_list)[i]
            if( is.logical(cc_set) ) {
                is_cell_control <- cc_set
                cc_set <- which(cc_set)
            } else {
                is_cell_control <- rep(FALSE, ncol(object))    
            }
            if(is.character(cc_set))
                cc_set <- which(cellNames(object) %in% cc_set)
            is_cell_control[cc_set] <- TRUE
            ## Construct data.frame for pData from this gene control set
            is_cell_control <- as.data.frame(is_cell_control)
            colnames(is_cell_control) <- paste0("is_cell_control_", set_name)
            if( exists("cell_controls_pdata") ) {
                cell_controls_pdata <- data.frame(cell_controls_pdata, 
                                                  is_cell_control)
            } else
                cell_controls_pdata <- is_cell_control
        }
    }
     
    ## Check column names and get cell controls across all sets
    if( n_sets_cell_controls == 1 ) {
        colnames(cell_controls_pdata) <- "is_cell_control"
    } else {
        ## Combine all cell controls
        is_cell_control <- apply(cell_controls_pdata, 1, any)
        cell_controls_pdata <- cbind(cell_controls_pdata, is_cell_control)
    }
    
    ## Add cell-level QC metrics to pData
    new_pdata <- as.data.frame(pData(object))
    ### Remove columns to be replaced
    to_replace <- colnames(new_pdata) %in% 
        c(colnames(gene_controls_pdata), colnames(cell_controls_pdata),
          colnames(log10_reads_from_gene_controls), "log10_reads_from_gene_controls")
    new_pdata <- new_pdata[, !to_replace]
    ### Add new QC metrics
    new_pdata$depth <- depth
    new_pdata$log10_depth <- log10(depth)
    new_pdata$coverage <- coverage
    new_pdata$filter_on_depth <- filter_on_depth
    new_pdata$filter_on_coverage <- filter_on_coverage
    new_pdata <- cbind(new_pdata, gene_controls_pdata, 
                       log10_reads_from_gene_controls)
    new_pdata$reads_from_biological_genes <- reads_from_biological_genes
    new_pdata$log10_reads_from_biological_genes <- 
        log10(reads_from_biological_genes+1) 
    new_pdata <- cbind(new_pdata, cell_controls_pdata)
    pData(object) <-  new("AnnotatedDataFrame", new_pdata)
    
    ## Add gene-level QC metrics to fData
    total_reads <- sum(counts(object))
    new_fdata <- as.data.frame(fData(object))
    ### Remove columns that are to be replaced
    to_replace <- colnames(new_fdata) %in% colnames(gene_controls_fdata)
    new_fdata <- new_fdata[, !to_replace]
    ### Add new QC information
    new_fdata$mean_exprs <- rowMeans(exprs(object))
    new_fdata$exprs_rank <- rank(rowMeans(exprs(object)))
    new_fdata$total_gene_reads <- rowSums(counts(object))
    new_fdata$log10_total_gene_reads <- log10(new_fdata$total_gene_reads+1)
    new_fdata$prop_total_reads <- rowSums(counts(object)) / total_reads
    new_fdata$n_cells_exprs <- rowSums(is_exprs(object))
    new_fdata <- cbind(new_fdata, gene_controls_fdata)
    fData(object) <- new("AnnotatedDataFrame", new_fdata)
    
    ## Ensure sample names are correct and return object
    sampleNames(object) <- colnames(exprs(object))
    object
}


################################################################################

#' Find most important principal components for a given variable
#'
#' @param object an SCESet object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param variable character scalar providing a variable name (column from 
#' \code{pData(object)}) for which to determine the most important PCs.
#'
#' @details Plot the top 5 most important PCs for a given variable. Importance 
#' here is defined as the average silhouette width (computed with 
#' \code{cluster::silhouette}) when cells are clustered by the values of the 
#' variable. If the variable is discrete, unique values of the variable define 
#' the clusters. If the variable is continuous, it is coerced into three 
#' discrete values, "low" (bottom quartile), "medium" (middle two quartiles) 
#' and "high" (upper quartile), and these levels are used to define clusters.
#' Cells are coloured by their cluster in the plot.
#'  
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' rownames(pd) <- pd$Cell
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset <- calculateQCMetrics(example_sceset)
#' findImportantPCs(example_sceset, variable="coverage")
#' 
findImportantPCs <- function(object, variable="coverage") {
    pca <- prcomp(t(exprs(object)), retx=TRUE, center=TRUE, scale.=TRUE)
    colnames(pca$x) <- paste("component", 1:ncol(pca$x))
    if(!(variable %in% colnames(pData(object))))
        stop("variable not found in pData(object). 
             Please make sure pData(object)[, variable] exists.")
    x <- pData(object)[, variable]
    x_na <- is.na(x)
    x <- x[!x_na]
    if(length(unique(x)) <= 1)
        stop("variable only has one unique value, so cannot determine important
             principal components.")
    ## Determine type of variable
    typeof_x <- .getTypeOfVariable(object, variable)
    if( typeof_x == "discrete" ) {
        ## If x is a discrete variable
        x_int <- as.integer(factor(x))
        ## Define colours for plotting
        colour_pal <- c("#1F77B499", "#FF7F0E99", "#2CA02C99", "#D6272899", 
                        "#9467BD99", "#8C564B99", "#E377C299", "#7F7F7F99", 
                        "#BCBD2299", "#17BECF99")
        ## Compute ave silhouette widths
        ave_sil_width <- .calculateSilhouetteWidth(x_int, pca$x[!x_na,])
    } else {
        ## Define plotting colours
        colour_pal <- c("#F7FCF099", "#CCEBC599", "#7BCCC499", "#2B8CBE99",
                        "#08408199")
        ## If x is a continuous variable - turn into an ordinal variable with three
        ## values - low, medium and high - or five values based on quintiles
        x_int_1 <- cut(x, c(min(x)-1, quantile(x, probs=c(0.2, 0.4, 0.6, 0.8)), 
                                    max(x)+1), right=TRUE)
        x_int_1 <- as.integer(x_int_1)
        x_int_2 <- cut(x, c(min(x)-1, quantile(x, probs=c(0.25, 0.75)), 
                          max(x)+1), right=TRUE)
        x_int_2 <- as.integer(x_int_2)
        ## Compute ave silhouette widths
        ave_sil_width_1 <- .calculateSilhouetteWidth(x_int_1, pca$x[!x_na,])
        sum_top5_1 <- sum(sort(ave_sil_width_1, decreasing=TRUE)[1:5])
        ave_sil_width_2 <- .calculateSilhouetteWidth(x_int_2, pca$x[!x_na,])
        sum_top5_2 <- sum(sort(ave_sil_width_2, decreasing=TRUE)[1:5])
        ## see which number of clusters gives a larger sum for top 5 sil widths
        if( sum_top5_1 > sum_top5_2) {
            ave_sil_width <- ave_sil_width_1
            x_int <- x_int_1
        } else {
            ave_sil_width <- ave_sil_width_2
            colour_pal <- colour_pal[c(1, 3, 5)]
            x_int <- x_int_2
        }
    }
    ## Tidy up names and choose top 5 most important PCs for the variable
    names(ave_sil_width) <- colnames(pca$x)
    colnames(pca$x) <- paste0(colnames(pca$x), "\n(avg silh wdth ",
                              formatC(signif(ave_sil_width, digits=2), digits=2,
                                      format="fg", flag="#"), ")")
    top5 <- order(ave_sil_width, decreasing=TRUE)[1:5]
    ## Define colours for points
    Col <- rep("firebrick", nrow(pca$x))
    Col[!x_na] <- colour_pal[x_int]
    ## Plot these bad boys
    ## Get rid of any NA variable columns of PC so that top ordering aligns
    par(bty="n", col.lab="gray60")
    pairs(pca$x[, top5], pch=21, col="gray70", bg=Col)
    #     ave_sil_width
}


.calculateSilhouetteWidth <- function(x, mat) {
    ave_sil_width <- rep(NA, ncol(mat))
    for(i in 1:ncol(mat)) {
        si <- cluster::silhouette(x, dist(mat[,i]))
        ave_sil_width[i] <- summary(si)$avg.width
    }
    ave_sil_width
}


# Range of Silhouette Width    Interpretation
# 0.71-1.0	A strong structure has been found
# 0.51-0.70	A reasonable structure has been found
# 0.26-0.50	The structure is weak and could be artificial
# < 0.25	No substantial structure has been found

################################################################################

#' Plot the genes with the highest read counts
#'
#' @param object an SCESet object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param col_by_variable variable name (must be a column name of pData(object))
#' to be used to assign colours to cell-level values.
#' @param n numeric scalar giving the number of the most expressed features to 
#' show. Default value is 50.
#' @param drop_genes a character, logical or numeric vector indicating which 
#' genes (or features) to drop when producing the plot. For example, control 
#' genes might be dropped to focus attention on contribution from biological 
#' rather than synthetic genes.
#'
#' @details Plot the percentage of reads accounted for by the top n most highly
#' expressed genes across the dataset.
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' rownames(pd) <- pd$Cell
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset <- calculateQCMetrics(example_sceset)
#' plotHighestReadCounts(example_sceset, col_by_variable="coverage")
#' plotHighestReadCounts(example_sceset, col_by_variable="Mutation_Status")
#' 
plotHighestReadCounts <- function(object, col_by_variable="coverage", n=50,
                              drop_genes=NULL) {
    ## Check that variable to colour points exists
    if(!(col_by_variable %in% colnames(pData(object)))) {
        warning("col_by_variable not found in pData(object). 
             Please make sure pData(object)[, variable] exists. Colours will not be plotted.")
        plot_cols <- FALSE
    } else
        plot_cols <- TRUE
    x <- pData(object)[, col_by_variable]
    #     x_na <- is.na(x)
    #     x <- x[!x_na]
    ## Determine type of variable
    typeof_x <- .getTypeOfVariable(object, col_by_variable)
    ## Figure out which genes to drop
    if( !(is.null(drop_genes) | length(drop_genes) == 0) ) {
        if(is.character(drop_genes))
            drop_genes <- which(rownames(object) %in% drop_genes)
        if(is.logical(drop_genes))
            object <- object[!drop_genes,]
        else
            object <- object[-drop_genes,]
    }
    ## Compute QC metrics on the (possibly) modified SCESet object to make sure
    ## we have the relevant values for this set of genes
    object <- calculateQCMetrics(object)
    ## Find the most highly expressed genes in this dataset
    ### Order by total gene reads across whole dataset
    oo <- fData(object)$total_gene_reads %>% order(decreasing=TRUE) 
    fdata <- fData(object)
    fdata$Gene <- factor(rownames(object), levels = rownames(object)[rev(oo)])
    ## Determine percentage reads accounted for by top genes across all cells
    total_reads <- sum(counts(object))
    top50_pctage <- 100 * sum(fdata$total_gene_reads[oo[1:n]]) / total_reads
    ## Determine percentage of reads for top genes by cell
    df_pct_reads_by_cell <- (100 * t(counts(object)[oo[1:n],]) / 
                                 pData(object)$depth)
#     df_pct_reads_by_cell_2 <- 100 * apply(counts(object)[oo[1:n],], 1, 
#                                           function(x) x / pData(object)$depth)
#     
    ## Melt dataframe so it is conducive to ggplot
    df_pct_reads_by_cell_long <- reshape2::melt(df_pct_reads_by_cell)
    df_pct_reads_by_cell_long$Var2 <- factor(df_pct_reads_by_cell_long$Var2, 
                                             levels = rownames(object)[rev(oo[1:n])])
    ## Add colour variable information
    if(typeof_x=="discrete")
        df_pct_reads_by_cell_long$colour_by <- factor(x)
    else
        df_pct_reads_by_cell_long$colour_by <- x
    ## Make plot
    plot_most_expressed <- df_pct_reads_by_cell_long %>%
        ggplot(aes_string(y="Var2", x="value", colour="colour_by")) +
        geom_point(size = 180/n, alpha = 0.6, shape=124) +
        ggtitle(paste0("Top ", n, " genes account for ", 
                       format(top50_pctage, digits = 3), "% of all read counts")) +
        #         coord_flip() +
        ylab("Gene") +
        xlab("% of total read counts") +
        theme_bw() +
        theme(legend.position=c(1, 0), legend.justification=c(1, 0),
              axis.text.x=element_text(colour="gray35"), 
              axis.text.y=element_text(colour="gray35"),
              axis.title.x=element_text(colour="gray35"), 
              axis.title.y=element_text(colour="gray35"),
              title=element_text(colour="gray35"))
    ## Sort of colouring of points
    if(typeof_x=="discrete") {
        plot_most_expressed <- plot_most_expressed + 
            ggthemes::scale_colour_tableau(name=col_by_variable) +
            geom_point(aes_string(x="100*as.numeric(prop_total_reads)", y="Gene"), 
                       data=fdata[oo[1:n],], size = 200/n, colour="gray30", 
                       shape=21, fill="wheat")
    } else {
        plot_most_expressed <- plot_most_expressed + 
            scale_colour_gradient(name=col_by_variable, low="lightgoldenrod", 
                                  high="firebrick4", space="Lab") +
            geom_point(aes_string(x="100*as.numeric(prop_total_reads)", y="Gene"), 
                       data=fdata[oo[1:n],], size = 200/n, colour="gray30", 
                       shape=21, fill="aliceblue")
    }
    plot_most_expressed   
}


.getTypeOfVariable <- function(object, variable) {
    ## Extract variable
    x <- pData(object)[, variable]
    ## Get type
    if(is.character(x) | is.factor(x)) {
        typeof_x <- "discrete"
    } else {
        if(is.integer(x)) {
            if(length(unique(x)) > 10)
                typeof_x <- "continuous"
            else
                typeof_x <- "discrete"
        } else {
            if(is.numeric(x))
                typeof_x <- "continuous"
            else {
                x <- as.character(x)
                typeof_x <- "discrete"
                warning("Unrecognised variable type. Variable being coerced to discrete.
                        Please make sure pData(object)[, variable] is a proper discrete or continuous variable")
            }        
        }
    }
    typeof_x
}

################################################################################


#' Plot explanatory variables ordered by percentage of phenotypic variance explained
#'
#' @param object an SCESet object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param variables vectory of variable names (must be columns of pData(object))
#' to be plotted.
#'
#' @details Function produces a pairs plot of the explanatory variables ordered
#' by the percentage of gene expression variance (as measured by R-squared in a 
#' marginal linear model) explained by variable. Median percentage R-squared is 
#' reported on the plot for each variable. Discrete variables are coerced to a 
#' factor and plotted as integers with jittering. Variables with only one unique 
#' value are quietly ignored. 
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' rownames(pd) <- pd$Cell
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset <- calculateQCMetrics(example_sceset)
#' vars <- names(pData(example_sceset))[c(2:3, 5:14)]
#' plotExplanatoryVariables(example_sceset, variables=vars)
#' 
plotExplanatoryVariables <- function(object, variables=c("coverage", "depth")) {
    ## Check that variables are defined
    variables_to_plot <- NULL
    for(var in variables) {
        if(!(var %in% colnames(pData(object)))) {
            warning(paste("variable", var, "not found in pData(object). 
             Please make sure pData(object)[, variable] exists. This variable will not be plotted."))
        } else {
            if(length(unique(pData(object)[, var])) <= 1) {
                warning(paste("variable", var, "only has one unique value, so R^2 is not meaningful.
This variable will not be plotted."))
            } else {
            variables_to_plot <- c(variables_to_plot, var)
            }
        }
    }
    ## Initialise matrix to score R^2 values for each gene for each variable
    rsquared_mat <- matrix(NA, nrow=nrow(object), ncol=length(variables_to_plot))
    val_to_plot_mat <- matrix(NA, nrow=ncol(object), ncol=length(variables_to_plot))
    colnames(rsquared_mat) <- colnames(val_to_plot_mat) <- variables_to_plot
    rownames(rsquared_mat) <- rownames(object)
    rownames(val_to_plot_mat) <- colnames(object)
    for(var in variables_to_plot) {
        x <- pData(object)[, var]
        #     x_na <- is.na(x)
        #     x <- x[!x_na]
        ## Determine type of variable
        typeof_x <- .getTypeOfVariable(object, var)
        if(typeof_x=="discrete") {
            x <- factor(x)
            val_to_plot_mat[, var] <- jitter(as.integer(x))
        } else {
            val_to_plot_mat[, var] <- x
        }
        design <- model.matrix(~x)
        rsquared_mat[, var] <- .getRSquared(exprs(object), design)
    }
    ## Get median R^2 for each variable, add to labels and order by median R^2
    median_rsquared <- apply(rsquared_mat, 2, median)
    colnames(val_to_plot_mat) <- paste0(colnames(val_to_plot_mat), 
                                        "\n(med var expl = ", 
                                        formatC(signif(100*median_rsquared, digits=3), 
                                                digits=3, format="fg", flag="#"),
                                        "%)")
    oo_median <- order(median_rsquared, decreasing=TRUE)
    ## Plot these bad boys
    par(bty="n", col.lab="gray60")
    pairs(val_to_plot_mat[, oo_median], pch=21, col="gray60", bg="gray80",
          col.lab="red")
}


#' @importFrom limma lmFit
.getRSquared <- function(y, design) {
    fit <- limma::lmFit(y, design=design)
    sst <- rowSums(y^2)
    ssr <- sst - fit$df.residual * fit$sigma^2
    (ssr/sst)
}



################################################################################

#' Produce QC diagnostic plots
#'
#' @param object an SCESet object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param type character scalar providing type of QC plot to compute: 
#' "most-expressed" (showing genes with highest expression), "find-pcs" (showing
#' the most important principal components for a given variable), or 
#' "explanatory-variables" (showing a set of explanatory variables plotted 
#' against each other, ordered by marginal variance explained).
#' @param ... arguments passed to \code{plotHighestReadCounts}, 
#' \code{plotImportantPCs} and \code{plotExplanatoryVariables} as appropriate.
#'
#' @details Calculate useful quality control metrics to help with pre-processing
#' of data and identification of potentially problematic genes and cells.
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data=sc_example_cell_info)
#' rownames(pd) <- pd$Cell
#' example_sceset <- newSCESet(countData=sc_example_counts, phenoData=pd)
#' example_sceset <- calculateQCMetrics(example_sceset)
#' plotQC(example_sceset, type="most", col_by_variable="Mutation_Status")
#' plotQC(example_sceset, type="find", variable="coverage")
#' vars <- names(pData(example_sceset))[c(2:3, 5:14)]
#' plotQC(example_sceset, type="expl", variables=vars)
#' 
plotQC <- function(object, type="most-expressed", ...) {
    type <- match.arg(type, c("most-expressed", "find-pcs", 
                              "explanatory-variables"))
    if(type == "most-expressed") {
        plot_out <- plotHighestReadCounts(object, ...)
        print(plot_out)
        return(plot_out)
    }
    if(type == "find-pcs") {        
        findImportantPCs(object, ...)
    }
    if(type == "explanatory-variables") {
        plotExplanatoryVariables(object, ...)
    }
}

