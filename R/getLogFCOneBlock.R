#' Get logFC values for replicates within one experimental block
#'
#' @param data_matrix numeric matrix with expression values
#' (log-scale) for features (rows) and samples/cells (columns).
#' @param group character vector or factor defining the two groups to
#' be compared
#' @return matrix of logFC values for all pairwise comparisons of
#' samples in group1 and group2
#' @export
getLogFCOneBlock <- function(data_matrix, group) {
    ## Check that group has the correct number of levels: only two groups currently supported
    group <- as.factor(group)
    if(nlevels(group) > 2)
        stop("More than two groups to be compared. Can only compare between two groups, currently.")
    ## Define indices for groups 1 and 2
    idx1 <- which(group == levels(group)[1])
    idx2 <- which(group == levels(group)[2])
    ngroup1 <- length(idx1)
    ngroup2 <- length(idx2)
    ## Define number of genes, comparisons and FC matrix
    ngenes <- nrow(data_matrix)
    ncells <- ncol(data_matrix)
    k <- ngroup1 * ngroup2
    logFC <- matrix(NA, nrow = ngenes, ncol = k)
    ## Get logFC for all pairwise comparisons of genes in group1 and group2
    for(i in seq_len(ngroup1)) {
        idx_this_insertion <- ((i - 1) * ngroup2 + 1):(i * ngroup2)
        insert_this <- data_matrix[, idx1[i]] - data_matrix[, idx2]
        logFC[, idx_this_insertion] <- insert_this
    }
    ## Return logFC matrix
    logFC
}

