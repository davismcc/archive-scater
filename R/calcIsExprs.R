## A suite of functions for calculating necessary quantities for analysis

#' Calculate which gene are expressed in which cells using a count
#' threshold
#'
#' @param object an SCESet object with expression and/or count data
#' @param lowerDetectionLimit numeric scalar giving the minimum expression level
#' for an expression observation in a cell for it to qualify as expressed
#' @param exprs_data character scalar indicating whether the count data ("counts") or the 
#' transformed expression data ("exprs") should be used to define if an 
#' observation is expressed or not
#' @return a logical matrix indicating whether or not a gene in a particular 
#' cell is expressed.
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData=sc_example_counts)
#' is_exprs(example_sceset) <- calcIsExprs(example_sceset, lowerDetectionLimit=1,
#' exprs_data="exprs")
calcIsExprs <- function(object, lowerDetectionLimit = NULL, exprs_data = "counts")
{
    ## Check that args are appropriate
    exprs_data <- match.arg(exprs_data, c("counts", "exprs"))
    if ( exprs_data == "counts" ) {
        dat_matrix <- counts(object)
    }
    else {
        if ( exprs_data == "exprs" )
            dat_matrix <- exprs(object)
        else
            stop("Did not recognise 'exprs_data' argument. Should be either 'counts' or 'exprs'.")
    }
    ## Extract lowerDetectionLimit if not provided
    if( is.null(lowerDetectionLimit) )
        lowerDetectionLimit <- object@lowerDetectionLimit
    ## Decide which observations are above detection limit and return matrix
    isexprs <- dat_matrix > lowerDetectionLimit
    rownames(isexprs) <- rownames(dat_matrix)
    colnames(isexprs) <- colnames(dat_matrix)
    isexprs
}

