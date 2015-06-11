## Function to summarise expression across features


#' Summarise expression values across feature
#' 
#' Create a new \code{SCESet} with counts summarised at a different feature 
#' level. A typical use would be to summarise transcript-level counts at gene
#' level.
#' 
#' @param object an \code{SCESet} object.
#' @param use_as_exprs character string indicating which slot of the 
#' assayData from the \code{SCESet} object should be used as expression values. 
#' Valid options are \code{'exprs'} the expression slot, \code{'tpm'} the 
#' transcripts-per-million slot or \code{'fpkm'} the FPKM slot.
#' @param summarise_by character string giving the column of \code{fDat(object)}
#' that will be used as the features for which summarised expression levels are 
#' to be produced. Default is \code{'feature_id'}.
#'  \code{"exprs"}.
#' 
#' @details Only transcripts-per-million (TPM) and fragments per kilobase of 
#' exon per million reads mapped (FPKM) expression values should be aggregated 
#' across features. Since counts are not scaled by the length of the feature, 
#' expression in counts units are not comparable within a sample without 
#' adjusting for feature length. Thus, we cannot sum counts over a set of 
#' features to get the expression of that set (for example, we cannot sum counts
#' over transcripts to get accurate expression estimates for a gene). See the 
#' following link for a discussion of RNA-seq expression units by Harold Pimentel:
#' \url{https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/}.
#' 
#' @export
#' 
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' fd <- new("AnnotatedDataFrame", data = 
#' data.frame(gene_id = featureNames(example_sceset), 
#' feature_id = paste("feature", rep(1:500, each = 4), sep = "_")))
#' fData(example_sceset) <- fd
#' example_sceset_summarised <- 
#' summariseExprsAcrossFeatures(example_sceset, use_as_exprs = "counts")
#' example_sceset_summarised <- 
#' summariseExprsAcrossFeatures(example_sceset, use_as_exprs = "exprs")
#' 
summariseExprsAcrossFeatures <- function(object, use_as_exprs = "tpm", 
                                         summarise_by = "feature_id") {
    if ( !is(object, "SCESet") )
        stop("Object must be an SCESet")
    if ( !(summarise_by %in% colnames(fData(object))) )
        stop("The summarise_by argument is not a column of fData(object).")
    ## Define an expression matrix depending on which values we're using
    use_as_exprs <- match.arg(use_as_exprs, c("exprs", "tpm", "fpkm", "counts"))
    exprs_mat <- switch(use_as_exprs,
                        exprs = exprs(object),
                        tpm = tpm(object),
                        fpkm = fpkm(object),
                        counts = counts(object))
    if ( use_as_exprs == "exprs" && object@logged ) {
        exprs_mat <- 2 ^ exprs_mat - object@logExprsOffset
    }
    ## Use reshape2 to make a long version of the expression matrix
    tmp_exprs <- data.frame(feature = fData(object)[[summarise_by]], exprs_mat)
    tmp_exprs_long <- suppressMessages(reshape2::melt(tmp_exprs))
    exprs_new <- reshape2::acast(tmp_exprs_long, feature ~ variable, sum)
    cat("Collapsing expression to", nrow(exprs_new), "features.")
    ## Create a new SCESet object
    pd <- new("AnnotatedDataFrame", pData(object))
    fd <- new("AnnotatedDataFrame", data.frame(exprs_collapsed_to = rownames(exprs_new)))
    rownames(fd) <- rownames(exprs_new)
    if ( use_as_exprs == "exprs" && object@logged ) {
        exprs_new <- log2(exprs_new + object@logExprsOffset)
    }
    sce_out <- switch(use_as_exprs,
                      exprs = newSCESet(exprsData = exprs_new, phenoData = pd, 
                                      featureData = fd, logged = TRUE),
                      tpm = newSCESet(tpmData = exprs_new, phenoData = pd, 
                                    featureData = fd),
                      fpkm = newSCESet(fpkmData = exprs_new, phenoData = pd, 
                                     featureData = fd),
                      counts = newSCESet(countData = exprs_new, phenoData = pd, 
                                       featureData = fd))
    ## Summarise other data in the object if present
    if ( use_as_exprs != "counts" && !is.null(counts(object)) ) {
        tmp_exprs <- data.frame(feature = fData(object)[[summarise_by]], 
                                counts(object))
        tmp_exprs_long <- reshape2::melt(tmp_exprs)
        counts(sce_out) <- reshape2::acast(tmp_exprs_long, feature ~ variable, 
                                           sum)
    }
    if ( use_as_exprs != "tpm" && !is.null(tpm(object)) ) {
        tmp_exprs <- data.frame(feature = fData(object)[[summarise_by]], 
                                tpm(object))
        tmp_exprs_long <- reshape2::melt(tmp_exprs)
        tpm(sce_out) <- reshape2::acast(tmp_exprs_long, feature ~ variable, 
                                           sum)
    }
    if ( use_as_exprs != "fpkm" && !is.null(fpkm(object)) ) {
        tmp_exprs <- data.frame(feature = fData(object)[[summarise_by]], 
                                fpkm(object))
        tmp_exprs_long <- reshape2::melt(tmp_exprs)
        fpkm(sce_out) <- reshape2::acast(tmp_exprs_long, feature ~ variable, 
                                        sum)
    }
    if ( !is.null(cpm(object)) ) {
        tmp_exprs <- data.frame(feature = fData(object)[[summarise_by]], 
                                cpm(object))
        tmp_exprs_long <- reshape2::melt(tmp_exprs)
        cpm(sce_out) <- reshape2::acast(tmp_exprs_long, feature ~ variable, 
                                         sum)
    }
    ## Use gene symbols for rownames
    sce_out
}

