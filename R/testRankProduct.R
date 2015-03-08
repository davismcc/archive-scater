# Wrapper functions for the computation of p-values for rank-product tests
# Davis McCarthy, July 2014

#' Conduct rank-product test
#'
#' @param dfilt DGEList object containing data and so on that has been
#' filtered and expressed cells determined
#' @param nmax.genes integer giving the maximum number of genes for
#' which to compute exact rank-product p-values
#' @param max.rprod.exact integer giving the maximum rank-product for
#' which to compute exact p-values. Exact p-values are very slow to
#' compute for large rank-products, and converge to approximate
#' p-values.
#' @return data frame with mean rank, rank product, number of ways of
#' obtaining given rank product, approximate p-value and exact p-value
#' @export
testRankProduct <- function(dfilt, nmax.genes = 500, max.rprod.exact = 100000) {
    
    #### NEED TO ADJUST TO GENERIC FOR SCESet Class Object #####
    
#     navec <- rep(NA, nrow(dfilt))
#     rp.results <- data.frame(Ave.CPM = navec, Ave.logFC = navec, N = navec)
#     rownames(rp.results) <- rownames(dfilt)
#     rep.names <- unique(dfilt$cells.info$culture)
#     logFC.df <- matrix(NA, nrow = nrow(dfilt), ncol = 2*length(rep.names) + 1)
#     colnames(logFC.df) <- c(paste0("logFC.", unique(dfilt$cells.info$culture)), "Ave.logFC")
#     treatment.names <- unique(dfilt$cells.info$perturbed)
#     for(i in 1:nrow(rp.results)) {
#         ## Code for one gene to do get logFC in each replicate experiment
#         dofit <- TRUE
# 
#         isexpr <- dfilt$isexpr_cpm[i,]
#         culture <- dfilt$cells.info$culture[isexpr]
#         treatment <- dfilt$cells.info$perturbed[isexpr]
#         if( (length(unique(culture)) > 1) & (length(unique(treatment)) > 1) ) {
#             design <- model.matrix(~ culture + treatment)
#             colnames(design) <- c("(Intercept)", "Expt", "Metformin")
#         }
#         else if( (length(unique(culture)) > 1) & (length(unique(treatment)) <= 1) ) {
#             warning("Gene expressed in only one treatment. Cannot do test on treatment.")
#             rp.results$N[i] <- sum(isexpr)
#             dofit <- FALSE
#         } else if( length(unique(treatment)) > 1 ) {
#             warning(paste0("Gene ", rownames(rp.results)[i], " expressed in only one culture/experiment (culture ", unique(culture), "). Results should be treated with caution."))
#             design <- model.matrix(~ treatment)
#             colnames(design) <- c("(Intercept)", "Metformin")
#         } else {
#             warning("Gene expressed in only one culture/experiment and one treatment. No DE results possible.")
#             rp.results$N[i] <- sum(isexpr)
#             dofit <- FALSE
#         }
#         if( dofit ) { ## Test to see if we should actually do the testing or not
#             x <- dfilt$log2cpm.pc0.25[i, isexpr]
#             rp.results$Ave.CPM[i] <- mean(x)
#             for( j in 1:length(rep.names) ) {
#                 ave.expr.1 <- mean(x[perturbed == treatment.names[1]] & culture == rep.names[j])
#                 ave.expr.2 <- mean(x[perturbed == treatment.names[2]] & culture == rep.names[j])
#                 logFC.df[, j] <- NA ### FIX ME
#             }
#         }
#     }
## Compute average logFC and logFC-rankings for all genes
}

#' Conduct RankProduct test for a data matrix
#'
#' @param data_matrix numeric matrix with expression values
#' (log-scale) for features (rows) and samples/cells (columns).
#' @param group character vector or factor defining the two groups to
#' be compared
#' @param block character vector or factor defining the experimental
#' blocks within which the two groups are to be compared
#' @param delta character string indication approximate p-value to be
#' returned for rank-product results: "upper", "lower", or "geometric"
#' (default; gives approximate p-value)
#' @param compute_pvals logical, should p-values be computed? Default
#' is \code{FALSE} as for experiments with many genes or cells,
#' p-values will be extremely expensice or impossible to compute
#' @return results object summarising RankProduct results
#' @export
testRankProduct.matrix <- function(data_matrix, group, block = NULL, delta = "geometric", compute_pvals = FALSE) {
    ngenes <- nrow(data_matrix)
    group <- as.factor(group)
    block <- as.factor(block)
    ## If no blocking, then assume one source for the data --> easy case
    if(is.null(block) | nlevels(block) == 1) {
        message("Treating study as having one experimental block.")
        ## get logFCs
        logFC <- getLogFCOneBlock(data_matrix, group)
        ## Compute median logFC to report for each gene
        median_logFC <- Biobase::rowMedians(logFC)
        ## Get ranks for each pairwise comparison of cells across groups
        ## Take log-ranks: true RP is likely to cause numeric overflow
        lranks <- getRanksOneBlock(logFC, group, logged = TRUE)
        lranks_up_grp1 <- lranks[[1]]
        lranks_up_grp2 <- lranks[[2]]
    }
    else {
        message("Computing rank-products across ", nlevels(block), " experimental blocks.")
        blk_idx <- 1
        for(blk in levels(block)) {
            message("Block ", blk_idx, "...")
            select_blk <- block == blk
            this_grp <- group[select_blk]
            ## get logFCs
            this_logFC <- getLogFCOneBlock(data_matrix[, select_blk], this_grp)
            ## Get ranks for each pairwise comparison of cells across groups
            ## Take log-ranks: true RP is likely to cause numeric overflow
            this_lranks <- getRanksOneBlock(this_logFC, this_grp, logged = TRUE)
            if( blk_idx == 1 ) {
                logFC <- this_logFC
                lranks_up_grp1 <- this_lranks[[1]]
                lranks_up_grp2 <- this_lranks[[2]]
            }
            else {
                logFC <- cbind(logFC, this_logFC)
                lranks_up_grp1 <- cbind(lranks_up_grp1, this_lranks[[1]])
                lranks_up_grp2 <- cbind(lranks_up_grp2, this_lranks[[2]])
            }
            blk_idx <- blk_idx + 1
        }
    }
    ## Compute median logFC to report for each gene
    median_logFC <- Biobase::rowMedians(logFC)
    ## Get mean log(rank-products) for up-regulation in group1 and up-regulation in group2
    ## This is equivalent to the log of the geometric mean of the ranks
    lrp_up_grp1 <- rowMeans(lranks_up_grp1)
    lrp_up_grp2 <- rowMeans(lranks_up_grp2)
    ## Report the smallest mean log-rank-product for each gene
    lrp <- pmin(lrp_up_grp1, lrp_up_grp2)
    ## Report whether expression higher in group1 or group2
    label1 <- paste0(levels(group), collapse = ">")
    label2 <- paste0(rev(levels(group)), collapse = ">")
    direction <- rep(label1, length(lrp))
    direction[lrp_up_grp1 > lrp_up_grp2] <- label2
    ## Compute approximate p-values for geometric mean of rank-products, as if k = 1
    if(compute_pvals) {
        ## Compute the number of "replicates" if computing p-values
        nreps <- ncol(logFC)
        pvals <- calcRankProdPValsBounds(exp(lrp), ngenes, nreps, delta = delta)
    }
    else
        pvals <- rep(NA, ngenes)
    ## Define output
    out_table <- data.frame(Median_logFC = median_logFC, Dir_of_Effect = direction,
                            Geom_Mean_Rank = exp(lrp), P_Value = pvals)
    rownames(out_table) <- rownames(data_matrix)
    out_table
}


