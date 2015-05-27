## Function to get feature annotation information from biomaRt


#' Get feature annotation information from Biomart
#' 
#' Use the \code{biomaRt} package to add feature annotation information to an 
#' \code{SCESet}. 
#' 
#' @param object an \code{SCESet} object
#' @param filters character vector defining the "filters" terms to pass to the
#' biomaRt::getBM function.
#' @param attributes character vector defining the biomaRt attributes to pass to
#' the \code{attributes} argument of \code{\link[biomaRt]{getBM}}.
#' @param feature_symbol character string defining the biomaRt attribute to be 
#' used to define the symbol to be used for each feature (which appears as the 
#' \code{feature_symbol} in fData(object), subsequently). Default is 
#' \code{"mgi_symbol"}, gene symbols for mouse. This should be changed if the 
#' organism is not Mus musculus!
#' @param feature_id character string defining the biomaRt attribute to be used 
#' to define the ID to be used for each feature (which appears as the 
#' \code{feature_id} in fData(object), subsequently). Default is 
#' \code{"ensembl_gene_id"}, Ensembl gene IDs for mouse. This should be changed 
#' if the organism is not Mus musculus!
#' @param biomart character string defining the biomaRt to be used. Default is 
#' \code{"ensembl"}.
#' @param dataset character string defining the biomaRt dataset to use. Default
#' is \code{"mmusculus_gene_ensembl"}, which should be changed if the organism 
#' is not the mouse!
#' @param host optional character string argument which can be used to select a
#' particular \code{"host"} from biomaRt to use. Useful for accessing archived
#' versions of biomaRt data. Default is \code{NULL}, in which case the current 
#' version of the biomaRt is used.
#' 
#' @details See the documentation for the biomaRt package, specifically for the
#' functions \code{useMart} and \code{getBM}, for information on what are 
#' permitted values for the filters, attributes, biomart, dataset and host 
#' arguments.
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' object <- getBMFeatureAnnos(object)
#' }
#' 
getBMFeatureAnnos <- function(object, filters="ensembl_transcript_id", 
                              attributes=c("ensembl_transcript_id", 
                                           "ensembl_gene_id", "mgi_symbol", 
                                           "chromosome_name", "transcript_biotype",
                                           "transcript_start", "transcript_end", 
                                           "transcript_count"), 
                              feature_symbol="mgi_symbol",
                              feature_id="ensembl_gene_id",
                              biomart="ensembl", 
                              dataset="mmusculus_gene_ensembl",
                              host=NULL) {
    ## Define Biomart Mart to use
    if( is.null(host) )
        bmart <- biomaRt::useMart(biomart=biomart, dataset=dataset)
    else 
        bmart <- biomaRt::useMart(biomart=biomart, dataset=dataset, host=host) 
    ## Define feature IDs from SCESet object
    feature_ids <- featureNames(object)
    ## Get annotations from biomaRt
    feature_info <- biomaRt::getBM(attributes=attributes, 
                                   filters=filters, 
                                   values=feature_ids, mart=bmart)
    ## Match the feature ids to the filters ids used to get info from biomaRt
    mm <- match(feature_ids, feature_info[[filters]])
    feature_info_full <- feature_info[mm, ]
    rownames(feature_info_full) <- feature_ids
    ## Define gene symbol and gene id
    feature_info_full$feature_symbol <- feature_info_full[[feature_symbol]]
    feature_info_full$feature_id <- feature_info_full[[feature_id]]
    ## Use rownames for gene symbol if gene symbol is missing
    na_symbol <- (is.na(feature_info_full$feature_symbol) | 
                      feature_info_full$feature_symbol == "")
    feature_info_full$feature_symbol[na_symbol] <- 
        rownames(feature_info_full)[na_symbol]
    ## Use rownames from SCESet object (feature IDs) for feature_id if na
    feature_info_full$feature_id[is.na(feature_info_full$feature_id)] <-
        rownames(feature_info_full)[is.na(feature_info_full$feature_id)]
    ## Need to drop any duplicated columns that we want to replace
    old_fdata <- fData(object)
    keep_cols <- !(colnames(old_fdata) %in% 
                       c("feature_symbol", "feature_id", attributes))
    if( sum(keep_cols) > 0) {
        colnames_old_fdata <- colnames(old_fdata)
        old_fdata <- as.data.frame(old_fdata[, keep_cols])
        colnames(old_fdata) <- colnames_old_fdata[keep_cols]
        new_fdata <- cbind(old_fdata, feature_info_full)
    } else 
        new_fdata <- feature_info_full
    ## Add new feature annotations to SCESet object
    fData(object) <- new_fdata
    ## Return SCESet object
    object
}

