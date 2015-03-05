## A suite of functions for calculating necessary quantities for analysis

#' Calculate which gene are expressed in which cells using a count
#' threshold
#'
#' @param object a DGEList object with count data
#' @param min_gene_count_in_cell integer giving the minimum gene count
#' in a cell for it to qualify as expressed
#' @return a DGEList object with an 'isexpr_counts' object added, a
#' matrix of TRUE/FALSE values indicating whether or not a gene in a
#' particular cell is expressed.
#' @export
#' @examples
#' d <- edgeR::DGEList(counts = matrix(rpois(100, lambda = 20), nrow = 10))
#' calcIsExprCounts(d)
calcIsExprCounts <- function(object, min_gene_count_in_cell = 5) {
    isexpr <- object$counts >= min_gene_count_in_cell
    object$isexpr_counts <- isexpr
    return(object)
}


#' Calculate which gene are expressed in which cells using a
#' counts-per-million threshold
#'
#' @param object a DGEList object with count data
#' @param min_gene_count_in_cell integer giving the minimum gene count
#' in a cell for it to qualify as expressed
#' @param normalized.lib.sizes logical passed to cpm(); whether or not
#' to use normalized library sizes to compute counts-per-million
#' @return a DGEList object with an 'isexpr_cpm' object added, a
#' matrix of TRUE/FALSE values indicating whether or not a gene in a
#' particular cell is expressed.
#' @details Could incorporate normalization factors when computing
#' counts-per-million, but this is not implemented yet.
#' @export
#' @examples
#' d <- edgeR::DGEList(counts = matrix(rpois(100, lambda = 20), nrow = 10))
#' d$cells.info <- data.frame(treatment = rep(c("A", "B"), each = 5), 
#' lib_size = colSums(d$counts))
#' calcIsExprCPM(d)
calcIsExprCPM <- function(object, min_gene_count_in_cell = 5, 
                          normalized.lib.sizes = FALSE) {
    cpm_thresh <- min_gene_count_in_cell / object$samples$lib.size * 1e06
    if( is.null(object$cpm) )
        object$cpm <- edgeR::cpm(object, normalized.lib.sizes = 
                                     normalized.lib.sizes)
    isexpr <- t(t(object$cpm) > cpm_thresh)
    object$isexpr_cpm <- isexpr
    object$cells.info$cpm_thresh <- cpm_thresh
    return(object)
}

