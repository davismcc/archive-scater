## Convenience function for computing QC metrics and adding to pData & fData
### This file contains definitions for the following functions:
### * calculateQCMetrics
### * findImportantPCs
### * plotExplanatoryVariables
### * plotHighestExprs
### * plotQC
### 
### * .calculateSilhouetteWidth
### * .getRSquared
### * .getTypeOfVariable

#### Some ideas for adding to QC metrics - proportion of library from top x features
#     subset(seq_real_estate_long, feature == 10) %>% 
#         ggplot(aes(x = Proportion_Library, colour = culture)) + 
#         geom_density(aes(fill = culture), alpha = 0.5) + 
#         facet_wrap(~perturbed, ncol=2) + 
#         theme_igray(16) + scale_colour_tableau() + scale_fill_tableau() + 
#         xlab("Cumulative proportion of library from 10 most expressed features") + 
#         ylab("Density")

#     subset(seq_real_estate_long, feature == 50) %>% 
#         ggplot(aes(x = norm.lib.size, y = Proportion_Library, 
#                    colour = culture)) + geom_point(size = 4, alpha = 0.7) + geom_rug(alpha = 0.7) + 
#         facet_wrap(~perturbed, ncol = 2) + theme_igray(16) + scale_colour_tableau() + 
#         scale_fill_tableau() + ylab("Cumulative proportion of library from 50 most expressed features") + 
#         xlab("Normalised library size")

#' Calculate QC metrics
#'
#' @param object an SCESet object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param feature_controls a character vector of feature names, or a logical vector, 
#' or a numeric vector of indices used to identify feature controls (for example,
#' ERCC spike-in genes, mitochondrial genes, etc).
#' @param cell_controls a character vector of cell (sample) names, or a logical 
#' vector, or a numeric vector of indices used to identify cell controls (for 
#' example, blank wells or bulk controls).
#' @param nmads numeric scalar giving the number of median absolute deviations to be 
#' used to flag potentially problematic cells based on depth (total number of 
#' counts for the cell, or library size) and coverage (number of features with 
#' non-zero expression). Default value is 5.
#' @param pct_feature_controls_threshold numeric scalar giving a threshold for 
#' percentage of expression values accounted for by feature controls. Used as to 
#' flag cells that may be filtered based on high percentage of expression from 
#' feature controls.
#'
#' @details Calculate useful quality control metrics to help with pre-processing
#' of data and identification of potentially problematic features and cells.
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
calculateQCMetrics <- function(object, feature_controls=NULL, cell_controls=NULL, 
                               nmads=5, pct_feature_controls_threshold = 80) {
    ## We must have an SCESet object
    if (!is(object, "SCESet"))
        stop("object must be an SCESet object.")
    
    ## Compute cell-level metrics
    if ( is.null(is_exprs(object)) ) {
        if (is.null(counts(object))) {
            stop("Please define is_exprs(object). E.g. use is_exprs(object) <- exprs(object) > 0.1")
        } else {            
            warning("is_exprs(object) is null. Defining 'is_exprs' using count data and 
a lower count threshold of 0.")   
            isexprs <- calcIsExprs(object, lowerDetectionLimit = 0)
            rownames(isexprs) <- rownames(counts(object))
            colnames(isexprs) <- colnames(counts(object))
            is_exprs(object) <- isexprs
        }       
    }    
    
    ## See what versions of the expression data are available in the object
    exprs_mat <- exprs(object)
    if ( object@logged )
        exprs_mat <- 2 ^ exprs_mat - object@logExprsOffset
    counts_mat <- counts(object)
    tpm_mat <- tpm(object)
    fpkm_mat <- fpkm(object)
    
    ## Compute coverage and find outliers
    coverage <- colSums(is_exprs(object))
    mad_coverage <- mad(coverage)
    med_coverage <- median(coverage)
    keep_coverage <- c(med_coverage - 5*mad_coverage, 
                       med_coverage + 5*mad_coverage)
    filter_on_coverage <- (findInterval(coverage, keep_coverage) != 1)
    ## Compute depth if counts are present
    if ( !is.null(counts_mat) )
        depth <- colSums(counts_mat)
    else 
        depth <- colSums(exprs_mat)
    mad_depth <- mad(log10(depth))
    med_depth <- median(log10(depth))
    keep_depth <- c(med_depth - 5*mad_depth, med_depth + 5*mad_depth)
    filter_on_depth <- (findInterval(log10(depth), keep_depth) != 1)
    
    ## Contributions from control features
    ### Determine if vector or list
    if ( is.null(feature_controls) | length(feature_controls) == 0 ) {
        exprs_from_feature_controls <- pct_exprs_from_feature_controls <- 
            rep(0, ncol(object))
        feature_controls_pdata <- data.frame(exprs_from_feature_controls,
                                          pct_exprs_from_feature_controls) 
        if ( !is.null(counts_mat) ) {
            counts_from_feature_controls <- pct_counts_from_feature_controls <- 
                rep(0, ncol(object))
            feature_controls_pdata <- cbind(feature_controls_pdata, 
                                         data.frame(counts_from_feature_controls,
                                              pct_counts_from_feature_controls))
        }
        if ( !is.null(tpm_mat) ) {
            tpm_from_feature_controls <- pct_tpm_from_feature_controls <- 
                rep(0, ncol(object))
            feature_controls_pdata <- cbind(feature_controls_pdata,
                                         data.frame(tpm_from_feature_controls,
                                                    pct_tpm_from_feature_controls))
        }
        if ( !is.null(fpkm_mat) ) {
            fpkm_from_feature_controls <- pct_fpkm_from_feature_controls <- 
                rep(0, ncol(object))
            feature_controls_pdata <- cbind(feature_controls_pdata,
                                         data.frame(fpkm_from_feature_controls,
                                                    pct_fpkm_from_feature_controls))
        }
        is_feature_control <- rep(FALSE, nrow(object))    
        feature_controls_fdata <- data.frame(is_feature_control)
        n_sets_feature_controls <- 1
    } else {
        if ( is.list(feature_controls) ) {
            feature_controls_list <- feature_controls
            n_sets_feature_controls <- length(feature_controls)
        }
        else {
            feature_controls_list <- list(feature_controls)
            n_sets_feature_controls <- 1
        }
        ## Cycle through the feature_controls list and add QC info
        for (i in seq_len(length(feature_controls_list)) ) {
            gc_set <- feature_controls_list[[i]]
            set_name <- names(feature_controls_list)[i]
            if ( is.logical(gc_set) ) {
                is_feature_control <- gc_set
                gc_set <- which(gc_set)
            } else {
                is_feature_control <- rep(FALSE, nrow(object))    
            }
            if (is.character(gc_set))
                gc_set <- which(rownames(object) %in% gc_set)
            ## Get summaries from exprs for feature control set and construct df
#             exprs_from_feature_controls <- colSums(exprs_mat[gc_set,])
#             pct_exprs_from_feature_controls <- (100 * exprs_from_feature_controls / 
#                 colSums(exprs_mat))
#             filter_on_pct_exprs_from_feature_controls <- 
#                 (pct_exprs_from_feature_controls > pct_feature_controls_threshold)
#             df_pdata_this <- data.frame(exprs_from_feature_controls,
#                                         pct_exprs_from_feature_controls,
#                                         filter_on_pct_exprs_from_feature_controls)
            df_pdata_this <- .get_qc_metrics_from_exprs_mat(
                exprs_mat, gc_set, pct_feature_controls_threshold, 
                calc_top_features = (n_sets_feature_controls == 1), 
                exprs_type = "exprs")
            if ( !is.null(counts_mat) ) {
                df_pdata_counts <- .get_qc_metrics_from_exprs_mat(
                    counts_mat, gc_set, pct_feature_controls_threshold,
                    calc_top_features = (n_sets_feature_controls == 1), 
                    exprs_type = "counts")
                df_pdata_this <- cbind(df_pdata_this, df_pdata_counts)
            }            
            if ( !is.null(tpm_mat) ) {
                df_pdata_tpm <- .get_qc_metrics_from_exprs_mat(
                    tpm_mat, gc_set, pct_feature_controls_threshold,
                    calc_top_features = (n_sets_feature_controls == 1), 
                    exprs_type = "tpm")
                df_pdata_this <- cbind(df_pdata_this, df_pdata_tpm)
            }  
            if ( !is.null(fpkm_mat) ) {
                df_pdata_fpkm <- .get_qc_metrics_from_exprs_mat(
                    fpkm_mat, gc_set, pct_feature_controls_threshold,
                    calc_top_features = (n_sets_feature_controls == 1), exprs_type = "fpkm")
                df_pdata_this <- cbind(df_pdata_this, df_pdata_fpkm)
            } 
            is_feature_control[gc_set] <- TRUE
            ## Define number of feature controls expressed
            n_detected_feature_controls <- colSums(is_exprs(object)[gc_set,])
            df_pdata_this$n_detected_feature_controls <- n_detected_feature_controls
            colnames(df_pdata_this) <- paste(colnames(df_pdata_this), set_name, 
                                             sep = "_")
            if ( exists("feature_controls_pdata") )
                feature_controls_pdata <- cbind(feature_controls_pdata, df_pdata_this)
            else
                feature_controls_pdata <- df_pdata_this
            ## Construct data.frame for fData from this feature control set
            df_fdata_this <- data.frame(is_feature_control)
            colnames(df_fdata_this) <- paste(colnames(df_fdata_this), set_name, 
                                          sep = "_")
            if ( exists("feature_controls_fdata") )
                feature_controls_fdata <- cbind(feature_controls_fdata, df_fdata_this)
            else 
                feature_controls_fdata <- df_fdata_this
        }
    }
    
    ## Fix column names and define feature controls across all control sets
    if ( n_sets_feature_controls == 1 ) {
        colnames(feature_controls_pdata) <- gsub("_$", "", 
                                              colnames(feature_controls_pdata))
        colnames(feature_controls_fdata) <- gsub("_$", "", 
                                              colnames(feature_controls_fdata))
    } else {
        ## Combine all feature controls
        is_feature_control <- apply(feature_controls_fdata, 1, any)
        feature_controls_fdata <- cbind(feature_controls_fdata, is_feature_control)
        ## Compute metrics using all feature controls
        df_pdata_this <- .get_qc_metrics_from_exprs_mat(
            exprs_mat, is_feature_control, pct_feature_controls_threshold,
            calc_top_features = TRUE, exprs_type = "exprs")
#         exprs_from_feature_controls <- colSums(exprs_mat[is_feature_control,])
#         pct_exprs_from_feature_controls <- (100 * exprs_from_feature_controls / 
#                                              colSums(exprs_mat))
#         filter_on_pct_exprs_from_feature_controls <- 
#             (pct_exprs_from_feature_controls > pct_feature_controls_threshold)
#         df_pdata_this <- data.frame(exprs_from_feature_controls,
#                                     pct_exprs_from_feature_controls,
#                                     filter_on_pct_exprs_from_feature_controls)
        if ( !is.null(counts_mat) ) {
            df_pdata_counts <- .get_qc_metrics_from_exprs_mat(
                counts_mat, is_feature_control, pct_feature_controls_threshold,
                calc_top_features = TRUE, exprs_type = "counts")
            df_pdata_this <- cbind(df_pdata_this, df_pdata_counts)
        }            
        if ( !is.null(tpm_mat) ) {
            df_pdata_tpm <- .get_qc_metrics_from_exprs_mat(
                tpm_mat, is_feature_control, pct_feature_controls_threshold,
                calc_top_features = TRUE, exprs_type = "tpm")
            df_pdata_this <- cbind(df_pdata_this, df_pdata_tpm)
        }  
        if ( !is.null(fpkm_mat) ) {
            df_pdata_fpkm <- .get_qc_metrics_from_exprs_mat(
                fpkm_mat, is_feature_control, pct_feature_controls_threshold,
                calc_top_features = TRUE, exprs_type = "fpkm")
            df_pdata_this <- cbind(df_pdata_this, df_pdata_fpkm)
        } 
        n_detected_feature_controls <- colSums(is_exprs(object)[is_feature_control,])
        df_pdata_this$n_detected_feature_controls <- n_detected_feature_controls
        feature_controls_pdata <- cbind(feature_controls_pdata, df_pdata_this)
    }
    
    
    ## Define counts from endogenous features
    qc_pdata <- feature_controls_pdata
    qc_pdata$exprs_from_endogenous_features <- colSums(exprs_mat) - 
        feature_controls_pdata$exprs_from_feature_controls
    if ( !is.null(counts_mat) ) {
        qc_pdata$counts_from_endogenous_features <- depth - 
            feature_controls_pdata$counts_from_feature_controls
    }
    if ( !is.null(tpm_mat) ) {
        qc_pdata$tpm_from_endogenous_features <- colSums(tpm_mat) - 
            feature_controls_pdata$tpm_from_feature_controls
    }
    if ( !is.null(fpkm_mat) ) {
        qc_pdata$fpkm_from_endogenous_features <- colSums(fpkm_mat) - 
            feature_controls_pdata$fpkm_from_feature_controls
    }
    ## Define log10 read counts from feature controls
    cols_to_log <- grep("^counts_from_|^exprs_from_|^tpm_from_|^fpkm_from_", 
                        colnames(qc_pdata))
    log10_cols <- log10(qc_pdata[, cols_to_log, drop = FALSE] + 1)
    colnames(log10_cols) <- paste0("log10_", colnames(qc_pdata)[cols_to_log])
    ## Combine into a big pdata object
    qc_pdata <- cbind(qc_pdata, log10_cols)
    
    ## Define cell controls
    ### Determine if vector or list
    if ( is.null(cell_controls) | length(cell_controls) == 0 ) {
        is_cell_control <- rep(FALSE, ncol(object))     
        cell_controls_pdata <- data.frame(is_cell_control)
        n_sets_cell_controls <- 1
    } else {
        if ( is.list(cell_controls) ) {
            cell_controls_list <- cell_controls
            n_sets_cell_controls <- length(cell_controls)
        }
        else {
            cell_controls_list <- list(cell_controls)
            n_sets_cell_controls <- 1
        }
        for (i in seq_len(n_sets_cell_controls) ) {
            cc_set <- cell_controls_list[[i]]
            set_name <- names(cell_controls_list)[i]
            if ( is.logical(cc_set) ) {
                is_cell_control <- cc_set
                cc_set <- which(cc_set)
            } else {
                is_cell_control <- rep(FALSE, ncol(object))    
            }
            if (is.character(cc_set))
                cc_set <- which(cellNames(object) %in% cc_set)
            is_cell_control[cc_set] <- TRUE
            ## Construct data.frame for pData from this feature control set
            is_cell_control <- as.data.frame(is_cell_control)
            colnames(is_cell_control) <- paste0("is_cell_control_", set_name)
            if ( exists("cell_controls_pdata") ) {
                cell_controls_pdata <- data.frame(cell_controls_pdata, 
                                                  is_cell_control)
            } else
                cell_controls_pdata <- is_cell_control
        }
    }
     
    ## Check column names and get cell controls across all sets
    if ( n_sets_cell_controls == 1 ) {
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
        c(colnames(qc_pdata), colnames(cell_controls_pdata))
    new_pdata <- new_pdata[, !to_replace, drop = FALSE]
    ### Add new QC metrics
    new_pdata$depth <- depth
    new_pdata$log10_depth <- log10(depth)
    new_pdata$filter_on_depth <- filter_on_depth
    new_pdata$coverage <- coverage
    new_pdata$filter_on_coverage <- filter_on_coverage
    new_pdata$pct_dropout <- 100 * colSums(!is_exprs(object)) / nrow(object)
    new_pdata <- cbind(new_pdata, qc_pdata, cell_controls_pdata)
    pData(object) <-  new("AnnotatedDataFrame", new_pdata)
    
    ## Add feature-level QC metrics to fData
    new_fdata <- as.data.frame(fData(object))
    ### Remove columns that are to be replaced
    to_replace <- colnames(new_fdata) %in% colnames(feature_controls_fdata)
    new_fdata <- new_fdata[, !to_replace, drop = FALSE]
    ### Add new QC information
    new_fdata$mean_exprs <- rowMeans(exprs(object))
    new_fdata$exprs_rank <- rank(rowMeans(exprs(object)))
    new_fdata$n_cells_exprs <- rowSums(is_exprs(object))
    total_exprs <- sum(exprs_mat)
    new_fdata$total_feature_exprs <- rowSums(exprs_mat)
    new_fdata$log10_total_feature_exprs <- log10(new_fdata$total_feature_exprs + 1)
    new_fdata$pct_total_exprs <- 100 * rowSums(exprs_mat) / total_exprs
    new_fdata$pct_dropout <- 100 * rowSums(!is_exprs(object)) / ncol(object)
  
    if ( !is.null(counts_mat) ) {
        total_counts <- sum(counts_mat)
        new_fdata$total_feature_counts <- rowSums(counts_mat)
        new_fdata$log10_total_feature_counts <- log10(new_fdata$total_feature_counts + 1)
        new_fdata$pct_total_counts <- 100 * rowSums(counts_mat) / total_counts
    }
    if ( !is.null(tpm_mat) ) {
        total_tpm <- sum(tpm_mat)
        new_fdata$total_feature_tpm <- rowSums(tpm_mat)
        new_fdata$log10_total_feature_tpm <- log10(new_fdata$total_feature_tpm + 1)
        new_fdata$pct_total_tpm <- 100 * rowSums(tpm_mat) / total_tpm
    }
    if ( !is.null(fpkm_mat) ) {
        total_fpkm <- sum(fpkm_mat)
        new_fdata$total_feature_fpkm <- rowSums(fpkm_mat)
        new_fdata$log10_total_feature_fpkm <- log10(new_fdata$total_feature_fpkm + 1)
        new_fdata$pct_total_fpkm <- 100 * rowSums(fpkm_mat) / total_fpkm
    }
    ## Add new fdata to object
    new_fdata <- cbind(new_fdata, feature_controls_fdata)
    fData(object) <- new("AnnotatedDataFrame", new_fdata)
    
    ## Ensure sample names are correct and return object
    sampleNames(object) <- colnames(exprs(object))
    object
}

.get_qc_metrics_from_exprs_mat <- function(exprs_mat, is_feature_control, 
                                           pct_feature_controls_threshold,
                                           calc_top_features = FALSE,
                                           exprs_type = "exprs") {
    ## Get total expression from feature controls
    exprs_from_feature_controls <- colSums(exprs_mat[is_feature_control,])
    ## Get % expression from feature controls
    pct_exprs_from_feature_controls <- (100 * exprs_from_feature_controls / 
                                         colSums(exprs_mat))
    ## Indicate whether or not to filter on percentage from controls
    filter_on_pct_exprs_from_feature_controls <- 
        (pct_exprs_from_feature_controls > pct_feature_controls_threshold)
    ## Make a data frame
    df_pdata_this <- data.frame(exprs_from_feature_controls,
                                pct_exprs_from_feature_controls,
                                filter_on_pct_exprs_from_feature_controls)
    if (calc_top_features) { ## Do we wnat to calculate exprs accounted for by 
        ## top features?
        ## Determine percentage of counts for top features by cell
        pct_exprs_by_cell <- (100 * t(t(exprs_mat) / colSums(exprs_mat)))
        pct_exprs_by_cell_sorted <- apply(pct_exprs_by_cell, 2, sort, 
                                          decreasing = TRUE)
        pct_exprs_top <- matrixStats::colCumsums(pct_exprs_by_cell_sorted[1:500,])
        pct_exprs_top_out <- t(pct_exprs_top[c(50, 100, 200),])
        colnames(pct_exprs_top_out) <- paste0("pct_exprs_from_top_", 
                                              c(50, 100, 200), "_features")
        df_pdata_this <- cbind(df_pdata_this, pct_exprs_top_out)
    }
    colnames(df_pdata_this) <- gsub("exprs", exprs_type, colnames(df_pdata_this))
    df_pdata_this
}



################################################################################

#' Find most important principal components for a given variable
#'
#' @param object an SCESet object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param variable character scalar providing a variable name (column from 
#' \code{pData(object)}) for which to determine the most important PCs.
#' @param plot_type character string, indicating which type of plot to produce.
#' Default, \code{"pairs-pcs"} produces a pairs plot for the top 5 PCs based on
#' their R-squared with the variable of interest. A value of 
#' \code{"pcs-vs-vars"} produces plots of the top PCs against the variable of 
#' interest.
#' @param theme_size numeric scalar providing base font size for ggplot theme.
#'
#' @details Plot the top 5 most important PCs for a given variable. Importance 
#' here is defined as the R-squared value from a linear model regressing each PC
#' onto the variable of interest. If the variable is discrete, unique values of 
#' the variable define  the clusters. If the variable is continuous, it is coerced into three 
#' discrete values, "low" (bottom quartile), "medium" (middle two quartiles) 
#' and "high" (upper quartile), and these levels are used as a factor in the 
#' linear model. Cells are coloured by their cluster in the plot.
#'  
#' @import viridis
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
findImportantPCs <- function(object, variable="coverage", 
                             plot_type = "pcs-vs-vars", theme_size = 10) {
    pca <- prcomp(t(exprs(object)), retx = TRUE, center = TRUE, scale. = TRUE)
    colnames(pca$x) <- paste("component", 1:ncol(pca$x))
    if (!(variable %in% colnames(pData(object))))
        stop("variable not found in pData(object). 
             Please make sure pData(object)[, variable] exists.")
    x <- pData(object)[, variable]
    x_na <- is.na(x)
    x <- x[!x_na]
    if (length(unique(x)) <= 1)
        stop("variable only has one unique value, so cannot determine important
             principal components.")
    ## Determine type of variable
    typeof_x <- .getTypeOfVariable(object, variable)
    if ( typeof_x == "discrete" ) {
        ## If x is a discrete variable
        x_int <- as.factor(x)
        ## Compute R-squared for each PC
        design <- model.matrix(~x_int)
    } else {
        ## If x is a continuous variable - turn into an ordinal variable with three
        ## values - low, medium and high - or five values based on quintiles
        x_int <- cut(x, c(min(x) - 1, quantile(x, probs = c(0.2, 0.4, 0.6, 0.8)), 
                                    max(x) + 1), right = TRUE)
        design <- model.matrix(~x)
    }
    ## Get R-squared for each PC for the variable of interest
    pca_r_squared <- .getRSquared(t(pca$x[!x_na,]), design)
    ## Tidy up names and choose top 5 most important PCs for the variable
    # names(ave_sil_width) <- colnames(pca$x)
    names(pca_r_squared) <- colnames(pca$x)
    colnames(pca$x) <- paste0(colnames(pca$x), "\n(R-squared ",
                              formatC(signif(pca_r_squared, digits = 2), 
                                      digits = 2, format = "fg", flag = "#"), ")")
    top5 <- order(pca_r_squared, decreasing = TRUE)[1:5]
    if ( plot_type == "pairs-pcs" ) {
#         ## Define colours for points
#         Col <- rep("firebrick", nrow(pca$x))
#         Col[!x_na] <- colour_pal[x_int]
#         ## Plot these bad boys
#         ## Get rid of any NA variable columns of PC so that top ordering aligns
#         par(bty = "n", col.lab = "gray60")
#         pairs(pca$x[, top5], pch = 21, col = "gray70", bg = Col)
        colour_by <- pData(object)[, variable]
        ## Generate a larger data.frame for pairs plot
        df_to_expand <- pca$x[, top5]
#         colnames(df_to_expand) <- colnames(pca$x)[, top5]
#         rownames(df_to_expand) <- sampleNames(object)
        names(df_to_expand) <- colnames(df_to_expand)
        gg1 <- .makePairs(df_to_expand)
        ## new data frame 
        df_to_plot_big <- data.frame(gg1$all, colour_by)
        # colnames(df_to_plot_big)[-c(1:4)] <- get("variable")
        ## pairs plot
        plot_out <- ggplot(df_to_plot_big, aes_string(x = "x", y = "y")) + 
            geom_point(aes_string(fill = "colour_by"), colour = "gray40", 
                       shape = 21, alpha = 0.65) +
            facet_grid(xvar ~ yvar, scales = "free") + 
            stat_density(aes_string(x = "x", 
                                    y = "(..scaled.. * diff(range(x)) + min(x))"),
                         data = gg1$densities, position = "identity", 
                         colour = "grey20", geom = "line") +
            xlab("") + 
            ylab("") +
            theme_bw(theme_size)
        if ( typeof_x == "discrete" ) {
            plot_out <- plot_out + 
                ggthemes::scale_fill_tableau(name = get("variable"))
        } else {
            plot_out <- plot_out +
                viridis::scale_fill_viridis(name = get("variable"))
        }
        return(plot_out)
    } else {
        top6 <- order(pca_r_squared, decreasing = TRUE)[1:6]
        df_to_plot <- reshape2::melt(pca$x[, top6])
        xvar <- pData(object)[, variable]
        df_to_plot$xvar <- rep(xvar, 6)
        pcs_vars_plot <- ggplot(df_to_plot, aes_string(x = "xvar", y = "value"),
                                colour = "black") +
            facet_wrap(~ Var2, nrow = 3, scales = "free_y") +
            xlab(variable) +
            ylab("Principal component value") + 
            theme_bw(theme_size)
        if ( typeof_x == "discrete") {
            pcs_vars_plot <- pcs_vars_plot + 
                geom_violin(fill = "aliceblue", colour = "gray60", 
                            alpha = 0.6, scale = "width") +
                geom_boxplot(width = 0.25, outlier.size = 0) +
                geom_dotplot(fill = "gray10", alpha = 0.6, shape = 21, 
                             binaxis = 'y', stackdir = 'center', dotsize = 1)
#                 geom_point(fill = "gray10", alpha = 0.5, shape = 21,
#                            position = position_jitter(height = 0, width = 0.2))
        } else {
            pcs_vars_plot <- pcs_vars_plot + 
                geom_point(fill = "gray10", alpha = 0.6, shape = 21) +
                stat_smooth(aes(group = 1), method = "lm", alpha = 0.3)
        }
        return(pcs_vars_plot)
    }
}



#' @importFrom limma lmFit
.getRSquared <- function(y, design) {
    ## Mean-centre rows to get correct R-squared values with the limma formula below
    y0 <- t(scale(t(y), center = TRUE, scale = FALSE))
    ## Get linear model fit
    fit <- limma::lmFit(y0, design = design)
    ## Compute total sum of squares
    sst <- rowSums(y0 ^ 2)
    ## Compute residual sum of squares
    ssr <- sst - fit$df.residual * fit$sigma ^ 2
    (ssr/sst)
}

# 
# .calculateSilhouetteWidth <- function(x, mat) {
#     ave_sil_width <- rep(NA, ncol(mat))
#     for (i in 1:ncol(mat)) {
#         si <- cluster::silhouette(x, dist(mat[,i]))
#         ave_sil_width[i] <- summary(si)$avg.width
#     }
#     ave_sil_width
# }

# Range of Silhouette Width    Interpretation
# 0.71-1.0	A strong structure has been found
# 0.51-0.70	A reasonable structure has been found
# 0.26-0.50	The structure is weak and could be artificial
# < 0.25	No substantial structure has been found

################################################################################

#' Plot the features with the highest expression values
#'
#' @param object an SCESet object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param col_by_variable variable name (must be a column name of pData(object))
#' to be used to assign colours to cell-level values.
#' @param n numeric scalar giving the number of the most expressed features to 
#' show. Default value is 50.
#' @param drop_features a character, logical or numeric vector indicating which 
#' features (e.g. genes, transcripts) to drop when producing the plot. For 
#' example, control genes might be dropped to focus attention on contribution 
#' from endogenous rather than synthetic genes.
#' @param use_as_exprs which slot of the \code{assayData} in the \code{object} 
#' should be used to define expression? Valid options are "counts" (default), 
#' "tpm", "fpkm" and "exprs".
#'
#' @details Plot the percentage of counts accounted for by the top n most highly
#' expressed features across the dataset.
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' rownames(pd) <- pd$Cell
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset <- calculateQCMetrics(example_sceset, feature_controls = 1:500)
#' plotHighestExprs(example_sceset, col_by_variable="coverage")
#' plotHighestExprs(example_sceset, col_by_variable="Mutation_Status")
#' 
plotHighestExprs <- function(object, col_by_variable = "coverage", n = 50,
                              drop_features = NULL, use_as_exprs = "counts") {
    ## Check that variable to colour points exists
    if (!(col_by_variable %in% colnames(pData(object)))) {
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
    ## Figure out which features to drop
    if ( !(is.null(drop_features) | length(drop_features) == 0) ) {
        if (is.character(drop_features))
            drop_features <- which(rownames(object) %in% drop_features)
        if (is.logical(drop_features))
            object <- object[!drop_features,]
        else
            object <- object[-drop_features,]
    }
    ## Compute QC metrics on the (possibly) modified SCESet object to make sure
    ## we have the relevant values for this set of features
    if ( !is.null(fData(object)$is_feature_control) )
        object <- calculateQCMetrics(
            object, feature_controls = fData(object)$is_feature_control)
    else
        object <- calculateQCMetrics(object)
   
    ## Define expression values to be used
    use_as_exprs <- match.arg(use_as_exprs, c("exprs", "tpm", "cpm", "fpkm", "counts"))
    exprs_mat <- switch(use_as_exprs,
                        exprs = exprs(object),
                        tpm = tpm(object),
                        cpm = cpm(object),
                        fpkm = fpkm(object),
                        counts = counts(object))
    if ( is.null(exprs_mat) ) {
        warning(paste0("The object does not contain ", use_as_exprs, " expression values. Using exprs(object) values instead."))
        exprs_mat <- exprs(object)
        use_as_exprs <- "exprs"
    }
    if ( use_as_exprs == "exprs" && object@logged )
        exprs_mat <- 2 ^ (exprs_mat) - object@logExprsOffset
    
    ## Find the most highly expressed features in this dataset
    ### Order by total feature counts across whole dataset
    fdata <- fData(object)
    oo <- order(fdata[[paste0("total_feature_", use_as_exprs)]], 
                decreasing = TRUE)
    fdata$feature <- factor(featureNames(object), levels = rownames(object)[rev(oo)])
    ## Check if is_feature_control is defined
    if ( is.null(fdata$is_feature_control) )
        fdata$is_feature_control <- rep(FALSE, nrow(fdata))
    if ( is.null(fdata$Feature) )
        fdata$Feature <- featureNames(object)
    
    ## Determine percentage expression accounted for by top features across all cells
    total_exprs <- sum(exprs_mat)
    total_feature_exprs <- fdata[[paste0("total_feature_", use_as_exprs)]]
    top50_pctage <- 100 * sum(total_feature_exprs[oo[1:n]]) / total_exprs
    ## Determine percentage of counts for top features by cell
    df_pct_exprs_by_cell <- (100 * t(exprs_mat[oo[1:n],]) / colSums(exprs_mat))
   
    ## Melt dataframe so it is conducive to ggplot
    df_pct_exprs_by_cell_long <- reshape2::melt(df_pct_exprs_by_cell)
    df_pct_exprs_by_cell_long$Var2 <- factor(df_pct_exprs_by_cell_long$Var2, 
                                             levels = rownames(object)[rev(oo[1:n])])
    ## Add colour variable information
    if (typeof_x == "discrete")
        df_pct_exprs_by_cell_long$colour_by <- factor(x)
    else
        df_pct_exprs_by_cell_long$colour_by <- x
    ## Make plot
    plot_most_expressed <- ggplot(df_pct_exprs_by_cell_long, 
                                  aes_string(y = "Var2", x = "value", 
                                             colour = "colour_by")) +
        geom_point(alpha = 0.6, shape = 124) +
        ggtitle(paste0("Top ", n, " account for ", 
                       format(top50_pctage, digits = 3), "% of total")) +
        ylab("Feature") +
        xlab(paste0("% of total ", use_as_exprs)) +
        theme_bw(8) +
        theme(legend.position = c(1, 0), legend.justification = c(1, 0),
              axis.text.x = element_text(colour = "gray35"), 
              axis.text.y = element_text(colour = "gray35"),
              axis.title.x = element_text(colour = "gray35"), 
              axis.title.y = element_text(colour = "gray35"),
              title = element_text(colour = "gray35"))
    ## Sort of colouring of points
    if (typeof_x == "discrete") {
        plot_most_expressed <- plot_most_expressed + 
            ggthemes::scale_colour_tableau(name = col_by_variable)
    } else {
        plot_most_expressed <- plot_most_expressed + 
            scale_colour_gradient(name = col_by_variable, low = "lightgoldenrod", 
                                  high = "firebrick4", space = "Lab")
    }
    plot_most_expressed + geom_point(
        aes_string(x = paste0("as.numeric(pct_total_", use_as_exprs, ")"), 
                   y = "Feature", fill = "is_feature_control"), 
        data = fdata[oo[1:n],], colour = "gray30", shape = 21) +
        scale_fill_manual(values = c("aliceblue", "wheat")) +
        guides(fill = guide_legend(title = "Feature control?"))
}


.getTypeOfVariable <- function(object, variable) {
    ## Extract variable
    x <- pData(object)[, variable]
    ## Get type
    if (is.character(x) || is.factor(x) || is.logical(x)) {
        typeof_x <- "discrete"
    } else {
        if (is.integer(x)) {
            if (length(unique(x)) > 10)
                typeof_x <- "continuous"
            else
                typeof_x <- "discrete"
        } else {
            if (is.numeric(x))
                typeof_x <- "continuous"
            else {
                x <- as.character(x)
                typeof_x <- "discrete"
                warning(paste0("Unrecognised variable type for ", variable, 
". Variable being coerced to discrete. Please make sure pData(object)[, variable] is a proper discrete or continuous variable"))
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
#' @param method character scalar indicating the type of plot to produce. If 
#' "density", the function produces a density plot of R-squared values for each 
#' variable when fitted as the only explanatory variable in a linear model. If 
#' "pairs", then the function produces a pairs plot of the explanatory variables
#' ordered by the percentage of feature expression variance (as measured by 
#' R-squared in a marginal linear model) explained.
#' @param use_as_exprs which slot of the \code{assayData} in the \code{object} 
#' should be used to define expression? Valid options are "exprs" (default), 
#' "tpm", "fpkm", "cpm", and "counts".
#' @param nvars_to_plot integer, the number of variables to plot in the pairs 
#' plot. Default value is 10.
#' @param min_marginal_r2 numeric scalar giving the minimal value required for
#' median marginal R-squared for a variable to be plotted. Only variables with a
#' median marginal R-squared strictly larger than this value will be plotted.
#' @param variables optional character vector giving the variables to be plotted.
#' Default is \code{NULL}, in which case all variables in \code{pData(object)} 
#' are considered and the \code{nvars_to_plot} variables with the highest median
#' marginal R-squared are plotted.
#' @param return_object logical, should an \code{SCESet} object with median 
#' marginal R-squared values added to \code{varMetadata(object)} be returned?
#' @param theme_size numeric scalar giving font size to use for the plotting 
#' theme
#' @param ... parameters to be passed to \code{\link{pairs}}.
#'
#' @details If the \code{method} argument is "pairs", then the function produces
#' a pairs plot of the explanatory variables ordered by the percentage of 
#' feature expression variance (as measured by R-squared in a marginal linear 
#' model) explained by variable. Median percentage R-squared is reported on the 
#' plot for each variable. Discrete variables are coerced to a factor and 
#' plotted as integers with jittering. Variables with only one unique value are 
#' quietly ignored. 
#' 
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
plotExplanatoryVariables <- function(object, method = "density", 
                                     use_as_exprs = "exprs", nvars_to_plot = 10,
                                     min_marginal_r2 = 0, variables = NULL, 
                                     return_object = FALSE, theme_size = 10, 
                                     ...) {
    ## Check method argument
    method <- match.arg(method, c("density", "pairs"))
    ## Checking arguments for expression values
    use_as_exprs <- match.arg(use_as_exprs, 
                              choices = c("exprs", "norm_exprs", "counts", 
                                          "norm_counts", "tpm", "norm_tpm", 
                                          "fpkm", "norm_fpkm", "cpm", "norm_cpm"))
    exprs_mat <- switch(use_as_exprs,
                        exprs = exprs(object),
                        norm_exprs = norm_exprs(object),
                        tpm = tpm(object),
                        norm_tpm = norm_tpm(object),
                        cpm = cpm(object),
                        norm_cpm = norm_cpm(object),
                        fpkm = fpkm(object),
                        norm_fpkm = norm_fpkm(object),
                        norm_counts = norm_counts(object))
    if ( is.null(exprs_mat) ) {
        warning(paste0("The object does not contain ", use_as_exprs, " expression values. Using exprs(object) values instead."))
        exprs_mat <- exprs(object)
        use_as_exprs <- "exprs"
    }
    
    ## Check that variables are defined
    if ( is.null(variables) ) {
        variables_to_plot <- varLabels(object)
    } else {
        variables_to_plot <- NULL
        for (var in variables) {
            if ( !(var %in% colnames(pData(object))) ) {
                warning(paste("variable", var, "not found in pData(object). 
                     Please make sure pData(object)[, variable] exists. This variable will not be plotted."))
            } else {
                variables_to_plot <- c(variables_to_plot, var)
            }
        }
    }
    variables_all <- varLabels(object)
      
    ## Initialise matrix to store R^2 values for each feature for each variable
    rsquared_mat <- matrix(NA, nrow = nrow(object), 
                           ncol = length(variables_all))
    val_to_plot_mat <- matrix(NA, nrow = ncol(object), 
                              ncol = length(variables_all))
    colnames(rsquared_mat) <- colnames(val_to_plot_mat) <- variables_all
    rownames(rsquared_mat) <- rownames(object)
    rownames(val_to_plot_mat) <- colnames(object)
    
    ## Get R^2 values for each feature and each variable
    for (var in variables_all) {
        if ( var %in% variables_to_plot ) {
            if (length(unique(pData(object)[, var])) <= 1) {
                warning(paste("variable", var, "only has one unique value, so R^2 is not meaningful.
This variable will not be plotted."))
                rsquared_mat[, var] <- NA
            } else {
                x <- pData(object)[, var]
                #     x_na <- is.na(x)
                #     x <- x[!x_na]
                ## Determine type of variable
                typeof_x <- .getTypeOfVariable(object, var)
                if ( typeof_x == "discrete" ) {
                    x <- factor(x)
                    val_to_plot_mat[, var] <- jitter(as.integer(x))
                } else {
                    val_to_plot_mat[, var] <- x
                }
                design <- model.matrix(~x)
                rsquared_mat[, var] <- .getRSquared(exprs_mat, design)
#                 rsq_base <- apply(exprs_mat, 1, function(y) {
#                     lm.first <- lm(y ~ -1 + design); summary(lm.first)$r.squared})
#                 all(abs(rsq_base - rsquared_mat[, var]) < 0.000000000001)
            }
        } else {
            rsquared_mat[, var] <- NA
        }
    }
    
    ## Get median R^2 for each variable, add to labels and order by median R^2
    median_rsquared <- apply(rsquared_mat, 2, median)
    oo_median <- order(median_rsquared, decreasing = TRUE)
    nvars_to_plot <- min(sum(median_rsquared > min_marginal_r2, na.rm = TRUE), 
                         nvars_to_plot)
    
    if ( method == "pairs") {
        #     median_rsquared <- apply(rsquared_mat, 2, quantile, probs = 0.1, 
        #                              na.rm = TRUE)
        colnames(val_to_plot_mat) <- paste0(colnames(val_to_plot_mat), 
                                            "\n(Med. R-sq = ", 
                                            formatC(signif(100*median_rsquared, 
                                                       digits = 3),  digits = 3,
                                                format = "fg", flag = "#"), "%)")
        
        ## Plot these bad boys
        par(bty = "o", col.lab = "gray60")
        pairs(val_to_plot_mat[, oo_median[1:nvars_to_plot]], pch = 21, 
              col = "gray60", bg = "gray80", col.lab = "red", ...)
    } else {
        df_to_plot <- suppressMessages(reshape2::melt(
            rsquared_mat[, oo_median[1:nvars_to_plot]]))
        colnames(df_to_plot) <- c("Feature", "Expl_Var", "R_squared")
        df_to_plot$Pct_Var_Explained <- 100 * df_to_plot$R_squared
        df_to_plot$Expl_Var <- factor(
            df_to_plot$Expl_Var, levels = colnames(rsquared_mat)[oo_median[1:nvars_to_plot]])
        plot_out <- ggplot(df_to_plot, aes_string(x = "Pct_Var_Explained", 
                                                  colour = "Expl_Var")) +
            geom_line(stat = "density", alpha = 0.7, size = 2, trim = TRUE) +
            geom_vline(xintercept = 1, linetype = 2) +
            scale_x_log10(breaks = 10 ^ (-3:2), labels = c(0.001, 0.01, 0.1, 1, 10, 100)) +
            xlab(paste0("% variance explained (log10-scale)")) +
            ylab("") +
            coord_cartesian(xlim = c(10 ^ (-3), 100)) +
            ggthemes::scale_color_tableau(name = "", palette = "tableau20")            
        if ( library(cowplot, logical.return = TRUE) )
            plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
        else
            plot_out <- plot_out + theme_bw(theme_size)
    }
    
    if ( return_object ) {
        ## Return object so that marginal R^2 are added to varMetadata
        varMetadata(object) <- data.frame(
            labelDescription = paste("Median marginal R-squared =", 
                                     median_rsquared))
        fdata <- fData(object)
        rsq_out <- rsquared_mat[, oo_median[1:nvars_to_plot]]
        colnames(rsq_out) <- paste0("Rsq_", colnames(rsq_out))
        fdata_new <- new("AnnotatedDataFrame", cbind(fdata, rsq_out))
        fData(object) <- fdata_new
        if ( method == "density" )
            print(plot_out)
        return(object)
    } else {
        if ( method == "density" )
            return(plot_out)
    }
}


# ### Emergency test code
# rsq <- fData(sce_simmons_filt_strict)[, grep("^Rsq", names(fData(sce_simmons_filt_strict)))]
# head(rsq)
# colSums(rsq > 0.99)
# 
# rsq[rsq$Rsq_batch > 0.99,]
# 
# plotExpression(sce_simmons_filt_strict, rownames(rsq)[rsq$Rsq_batch > 0.99], x = "batch",
#                ncol = 6, colour_by = "tissue_type", use_as_exprs = "norm_exprs")
# 
# design <- model.matrix(~sce_simmons_filt_strict$batch)
# lm.first <- lm(norm_exprs(sce_simmons_filt_strict)["ERCC-00002",] ~ -1 + design)
# summary(lm.first)$r.squared
# summary(lm.first)$residuals
# plot(lm.first)
# lm2 <- lm(norm_exprs(sce_simmons_filt_strict)["ERCC-00002",] ~ sce_simmons_filt_strict$batch)
# summary(lm2)$r.squared
# 
# rsq_limma <- .getRSquared(norm_exprs(sce_simmons_filt_strict), design)
# 
# rsq_base2 <- apply(norm_exprs(sce_simmons_filt_strict), 1, function(y) {
#     lm.first <- lm(y ~ sce_simmons_filt_strict$batch); summary(lm.first)$r.squared})
# 
# all(abs(rsq_base2 - rsq_limma) < 0.00000000000001)
# 
# rsq_base <- apply(norm_exprs(sce_simmons_filt_strict), 1, function(y) {
#     lm.first <- lm(y ~ -1 + design); summary(lm.first)$r.squared})
# all(abs(rsq_base - rsq$Rsq_batch) < 0.000000000001)


################################################################################
### Plot expression frequency vs mean for feature controls

#' Plot frequency of expression against mean expression level
#' 
#' @param object an \code{SCESet} object.
#' @param feature_set character, numeric or logical vector indicating a set of
#' features to use for the PCA. If character, entries must all be in 
#' \code{featureNames(object)}. If numeric, values are taken to be indices for 
#' features. If logical, vector is used to index features and should have length
#' equal to \code{nrow(object)}. If \code{NULL}, then the function checks if 
#' feature controls are defined. If so, then only feature controls are plotted,
#' if not, then all features are plotted.
#' @param shape (optional) numeric scalar to define the plotting shape. 
#' @param alpha (optional) numeric scalar (in the interval 0 to 1) to define the
#'  alpha level (transparency) of plotted points.
#' @param ... further arguments passed to \code{\link{plotMetadata}} (should 
#' only be \code{size}, if anythin).
#' 
#' @export
#' @examples 
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data=sc_example_cell_info)
#' rownames(pd) <- pd$Cell
#' ex_sceset <- newSCESet(countData=sc_example_counts, phenoData=pd)
#' ex_sceset <- calculateQCMetrics(ex_sceset)
#' plotExprsFreqVsMean(ex_sceset)
#' 
plotExprsFreqVsMean <- function(object, feature_set = NULL, shape = 1,
                                alpha = 0.7, ...) {
    if ( !is(object, "SCESet") )
        stop("Object must be an SCESet")
    if ( is.null(fData(object)$n_cells_exprs) || 
         is.null(fData(object)$mean_exprs)) {
        stop("fData(object) does not have both 'n_cells_exprs' and 'mean_exprs' columns. Try running 'calculateQCMetrics' on this object, and then rerun this command.") 
    }
    if ( !is.null(feature_set) && feature_set != "feature_controls" && 
         typeof(feature_set) == "character" ) {
        if ( !(all(feature_set %in% featureNames(object))) )
            stop("when the argument 'feature_set' is of type character and not 'feature_controls', all features must be in featureNames(object)")
    }
    if ( is.null(feature_set) ) {
        feature_set <- 1:nrow(object)
        x_lab <- "Mean expression level (all features)"
    } else {
        if ( length(feature_set) == 1 && feature_set == "feature_controls" ) {
            feature_set <- fData(object)$is_feature_control
            x_lab <- "Mean expression level (feature controls)"
        } else
            x_lab <- "Mean expression level (supplied feature set)"
    }
    
    ## Plot this 
    if ( any(fData(object)$is_feature_control[feature_set]) && 
         !all(fData(object)$is_feature_control[feature_set]) ) {
        plot_out <- plotFeatureData(object[feature_set,], 
                                    aes_string(x = "mean_exprs", 
                                               y = "n_cells_exprs", 
                                               colour = "is_feature_control",
                                               shape = "is_feature_control"),
                                    alpha = alpha, ...) +
            scale_shape_manual(values = c(1, 17)) +
            ylab("Frequency of expression (number of cells)") +
            xlab(x_lab)
    } else {
        plot_out <- plotFeatureData(object[feature_set,], 
                                    aes_string(x = "mean_exprs", 
                                               y = "n_cells_exprs"), 
                                    alpha = alpha, shape = shape, ...) +
            ylab("Frequency of expression (number of cells)") +
            xlab(x_lab)
    }
    plot_out
}




################################################################################

#' Produce QC diagnostic plots
#'
#' @param object an SCESet object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param type character scalar providing type of QC plot to compute: 
#' "highest-expression" (showing features with highest expression), "find-pcs" (showing
#' the most important principal components for a given variable), 
#' "explanatory-variables" (showing a set of explanatory variables plotted 
#' against each other, ordered by marginal variance explained), or 
#' "exprs-mean-vs-freq" (plotting the mean expression levels against the 
#' frequency of expression for a set of features).
#' @param ... arguments passed to \code{plotHighestExprs}, 
#' \code{plotImportantPCs}, \code{plotExplanatoryVariables} and 
#' \code{plotExprsMeanVsFreq} as appropriate.
#'
#' @details Calculate useful quality control metrics to help with pre-processing
#' of data and identification of potentially problematic features and cells.
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data=sc_example_cell_info)
#' rownames(pd) <- pd$Cell
#' example_sceset <- newSCESet(countData=sc_example_counts, phenoData=pd)
#' example_sceset <- calculateQCMetrics(example_sceset)
#' plotQC(example_sceset, type="high", col_by_variable="Mutation_Status")
#' plotQC(example_sceset, type="find", variable="coverage")
#' vars <- names(pData(example_sceset))[c(2:3, 5:14)]
#' plotQC(example_sceset, type="expl", variables=vars)
#' 
plotQC <- function(object, type = "highest-expression", ...) {
    type <- match.arg(type, c("highest-expression", "find-pcs", 
                              "explanatory-variables", "exprs-freq-vs-mean"))
    if (type == "highest-expression") {
        plot_out <- plotHighestExprs(object, ...)
        return(plot_out)
    }
    if (type == "find-pcs") {        
        plot_out <- findImportantPCs(object, ...)
        if ( !is.null(plot_out) )
            return(plot_out)
    }
    if (type == "explanatory-variables") {
        plot_out <- plotExplanatoryVariables(object, ...)
        if ( !is.null(plot_out) )
            return(plot_out)
    }
    if (type == "exprs-freq-vs-mean") {
        plot_out <- plotExprsFreqVsMean(object, ...)
        if ( !is.null(plot_out) )
            return(plot_out)
    }
    
}

