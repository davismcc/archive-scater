## A suite of functions for calculating necessary quantities for analysis

#' Calculate which gene are expressed in which cells using a count
#' threshold
#'
#' @param object an SCESet object with expression and/or count data
#' @param lowerDetectionLimit numeric scalar giving the minimum expression level
#' for an expression observation in a cell for it to qualify as expressed
#' @param expr_data character scalar indicating whether the count data ("counts") or the 
#' transformed expression data ("exprs") should be used to define if an 
#' observation is expressed or not
#' @return a logical matrix indicating whether or not a gene in a particular 
#' cell is expressed.
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData = sc_example_counts)
#' isExpr(example_sceset) <- calcIsExpr(example_sceset, lowerDetectionLimit = 1,
#' expr_data = "exprs")
calcIsExpr <- function(object, lowerDetectionLimit = NULL, expr_data = "counts")
{
    ## Check that args are appropriate
    expr_data <- match.arg(expr_data, c("counts", "exprs"))
    if ( expr_data == "counts" ) {
        dat_matrix <- counts(object)
    }
    else {
        if ( expr_data == "exprs" )
            dat_matrix <- exprs(object)
        else
            stop("Did not recognise 'expr_data' argument. Should be either 'counts' or 'exprs'.")
    }
    ## Extract lowerDetectionLimit if not provided
    if( is.null(lowerDetectionLimit) )
        lowerDetectionLimit <- object@lowerDetectionLimit
    ## Decide which observations are above detection limit and return matrix
    isexpr <- dat_matrix > lowerDetectionLimit
    rownames(isexpr) <- rownames(dat_matrix)
    colnames(isexpr) <- colnames(dat_matrix)
    isexpr
}

