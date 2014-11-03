## Scale expression values to have average amplitude of 1 in
## expressing cells

#' Scale expression values to have average amplitude of 1 in
#' expressing cells
#'
#' @param expr_matrix a numeric matrix containing expression values
#' counts.
#' @param is_expr a logical matrix of same dimensions as
#' \code{expr_matrix} indicating which cells are expressed for each
#' gene.
#' @details This function is intended to process RNA-Seq or ChIP-Seq
#' data prior to modeling with RankProduct. Written by Davis McCarthy.
#' 
#' @return A matrix object.
#' @export
#'
scaleAmp <- function(expr_matrix, is_expr) {
    ## Get average amplitude in expressing cells for each gene
    ave_expr <- rep(NA, nrow(expr_matrix))
    for( i in seq_len(nrow(expr_matrix)) ) {
        expr_ind <- is_expr[i, ]
        expr_val <- expr_matrix[i, ]
        ave_expr[i] <- mean(expr_val[expr_ind])        
    }
    ## Return scaled expression values
    expr_matrix / ave_expr
}
