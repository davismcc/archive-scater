#' Get ranks from logFCs for one block
#'
#' Get ranks both for up-regulation in group1 and up-regulation in
#' group2
#'
#' @param logFC_matrix numeric matrix of logFCs computed as group1 -
#' group2
#' @param group character vector or factor defining the experimental
#' groups being compared
#' @param logged boolean: return raw ranks if \code{FALSE} or
#' log(ranks) if \code{TRUE}
#' @export
getRanksOneBlock <- function(logFC_matrix, group, logged = FALSE) {
    group <- as.factor(group)
    label1 <- paste0(levels(group), collapse = "-")
    label2 <- paste0(rev(levels(group)), collapse = "-")
    ## We assume that logFCs have been computed such that group1 = levels(group)[1] and logFC is group1 - group2
    ## Ranks in ascending order, so make logFC negative so that large positive logFCs end up at the top of the ranking
    rank_matrix_up_grp1 <- apply(-logFC_matrix, 2, rank)
    ## Rank by up-regulation in group2
    rank_matrix_up_grp2 <- apply(logFC_matrix, 2, rank)
    if(logged) {
        rank_matrix_up_grp1 <- log(rank_matrix_up_grp1)
        rank_matrix_up_grp2 <- log(rank_matrix_up_grp2)
    }
    ranks <- list(rank_matrix_up_grp1, rank_matrix_up_grp2)
    names(ranks) <- c(label1, label2)
    ranks
}

