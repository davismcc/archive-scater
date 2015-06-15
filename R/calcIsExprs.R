## A suite of functions for calculating necessary quantities for analysis

#' Calculate which features are expressed in which cells using a threshold on 
#' observed counts, transcripts-per-million, counts-per-million, FPKM, or 
#' defined expression levels.
#'
#' @param object an SCESet object with expression and/or count data.
#' @param lowerDetectionLimit numeric scalar giving the minimum expression level
#' for an expression observation in a cell for it to qualify as expressed.
#' @param exprs_data character scalar indicating whether the count data 
#' (\code{"counts"}), the transformed expression data (\code{"exprs"}), 
#' transcript-per-million (\code{"tpm"}), counts-per-million (\code{"cpm"}) or 
#' FPKM (\code{"fpkm"}) should be used to define if an observation is expressed 
#' or not.
#' @return a logical matrix indicating whether or not a feature in a particular 
#' cell is expressed.
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData=sc_example_counts)
#' is_exprs(example_sceset) <- calcIsExprs(example_sceset, lowerDetectionLimit = 1,
#' exprs_data = "exprs")
calcIsExprs <- function(object, lowerDetectionLimit = NULL, exprs_data = "counts")
{
    if( !is(object, "SCESet") )
        stop("Object must be an SCESet.")
    ## Check that args are appropriate
    exprs_data <- match.arg(exprs_data, c("counts", "exprs", "tpm", "cpm", "fpkm"))
    dat_matrix <- switch(exprs_data,
                         counts = counts(object),
                         exprs = exprs(object),
                         tpm = tpm(object),
                         cpm = cpm(object),
                         fpkm = fpkm(object))   
    if ( is.null(dat_matrix) )
        stop(paste0("Tried to use ", exprs_data, " as expression data, but ", 
                    exprs_data, "(object) is null."))
#     
#     if ( exprs_data == "counts" ) {
#         dat_matrix <- counts(object)
#     }
#     else {
#         if ( exprs_data == "exprs" )
#             dat_matrix <- exprs(object)
#         else
#             
#     }
    ## Extract lowerDetectionLimit if not provided
    if ( is.null(lowerDetectionLimit) )
        lowerDetectionLimit <- object@lowerDetectionLimit
    ## Decide which observations are above detection limit and return matrix
    isexprs <- dat_matrix > lowerDetectionLimit
    rownames(isexprs) <- rownames(dat_matrix)
    colnames(isexprs) <- colnames(dat_matrix)
    isexprs
}

