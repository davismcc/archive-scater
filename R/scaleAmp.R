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
    expr_vals <- expr_matrix
    expr_vals[!is_expr] <- NA
    ave_expr <- rowMeans(expr_vals, na.rm = TRUE)
    ## ave_expr <- rep(NA, nrow(expr_matrix))
    ## for( i in seq_len(nrow(expr_matrix)) ) {
    ##     expr_ind <- is_expr[i, ]
    ##     expr_val <- expr_matrix[i, ]
    ##     ave_expr[i] <- mean(expr_val[expr_ind], na.rm = TRUE)
    ## }
    ## Scale expression values for expressed observations
    out <- expr_matrix / ave_expr
    ## Leave expression values for non-expressed observations unscaled
    ## out[!is_expr] <- expr_matrix[!is_expr]
    out
}


## Scale expression values to have unit variance in expressing cells

#' Scale expression values to have unit variance in expressing cells
#'
#' @param expr_matrix a numeric matrix containing expression values
#' counts.
#' @param is_expr a logical matrix of same dimensions as
#' \code{expr_matrix} indicating which cells are expressed for each
#' gene.
#' @param vars optional vector of variance values to apply for the
#' scaling of the expression values. If \code{NULL} (default), then
#' variance is computed for each gene only from the cells with
#' non-zero expression for the gene.
#' @details This function is intended to process RNA-Seq or ChIP-Seq
#' data prior to modeling with RankProduct. Written by Davis McCarthy.
#'
#' @return A matrix object.
#' @export
#'
scaleVar <- function(expr_matrix, is_expr, vars = NULL) {
    if( is.null(vars) ) {
        ## Get variance of expression values in expressing cells for each gene
        ## Use rowVars() from the matrixStats package
        expr_vals <- expr_matrix
        expr_vals[!is_expr] <- NA
        var_expr <- matrixStats::rowVars(expr_vals, na.rm = TRUE)
        ## If any variances are zero, convert to 1 to keep observations as is
        var_expr[var_expr == 0] <- 1
        ## If any variances are NA, similarly replace the NA with 1
        var_expr[is.na(var_expr)] <- 1
    }
    else {
        var_expr <- vars
    }
    ## Scale expression values for expressed observations
    out <- expr_matrix / sqrt(var_expr)
    ## Leave expression values for non-expressed observations unscaled
    ## out[!is_expr] <- expr_matrix[!is_expr]
    out
}





