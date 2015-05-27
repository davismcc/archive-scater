## A set of functions for calculating and summarising expression values

#' Calculate transcripts-per-million (TPM)
#' 
#' Calculate transcripts-per-million (TPM) values for expression from counts for
#' a set of features.
#' 
#' @param object an \code{SCESet} object
#' @param effective_length vector of class \code{"numeric"} providing the 
#' effective length for each feature in the \code{SCESet} object
#' @param calc_from character string indicating whether to compute TPM from 
#' \code{"counts"} or \code{"fpkm"}. Default is to use \code{"counts"}, in which 
#' case the \code{effective_length} argument must be supplied.
#' 
#' @return Matrix of TPM values.
#' 
#' @examples
#' \dontrun{
#' tpm(object) <- calculateTPM(object, effective_length, calc_from="counts")
#' }
calculateTPM <- function(object, effective_length=NULL, calc_from="counts") {
    calc_from <- match.arg(calc_from, c("counts", "fpkm"), several.ok=FALSE)
    if( calc_from=="fpkm" )
        tpm_to_add <- .fpkmToTpm(fpkm(object))
    else {
        if( is.null(effective_length) )
            stop("effective_length argument is required if computing TPM from counts")
        tpm_to_add <- .countToTpm(counts(object), effective_length)
    }
    ## Return TPM values
    tpm_to_add
}


#' Calculate fragments per kilobase of exon per million reads mapped (FPKM)
#' 
#' Calculate fragments per kilobase of exon per million reads mapped (FPKM) 
#' values for expression from counts for a set of features.
#' 
#' @param object an \code{SCESet} object
#' @param effective_length vector of class \code{"numeric"} providing the 
#' effective length for each feature in the \code{SCESet} object
#' 
#' @return Matrix of FPKM values.
#' 
#' @examples
#' \dontrun{
#' fpkm(object) <- calculateFPKM(object, effective_length)
#' }
calculateFPKM <- function(object, effective_length) {
    fpkm_to_add <- .countToFpkm(counts(object), effective_length)
    ## Return matrix of FPKM values
    fpkm_to_add
}



.countToTpm <- function(counts, eff_len) {
    ## Expecting a count matrix of ngenes x ncells
    ## can't have any zero counts, so expect to apply offset
    rate <- log(counts) - log(eff_len)
    denom <- log(colSums(exp(rate)))
    exp( t(t(rate) - denom) + log(1e6) )
}

.countToFpkm <- function(counts, eff_len) {
    ## Expecting a count matrix of ngenes x ncells
    ## can't have any zero counts, so expect to apply offset
    N <- colSums(counts)
    logfpkm <- log(counts) + log(1e9) - log(eff_len)
    logfpkm <- t(t(logfpkm) - log(N))
    exp( logfpkm )
}

.fpkmToTpm <- function(fpkm) {
    ## Expecting an FPKM matrix of ngenes x ncells
    exp( t(t(log(fpkm)) - log(colSums(fpkm))) + log(1e6) )
}

.countToEffCounts <- function(counts, len, eff_len) {
    counts * (len / eff_len)
}






