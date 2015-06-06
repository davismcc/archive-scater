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
#' 
#' @details Only transcripts-per-million (TPM) and fragments per kilobase of 
#' exon per million reads mapped (FPKM) expression values should be aggregated 
#' across features. Since counts are not scaled by the length of the feature, 
#' expression in counts units are not comparable within a sample without 
#' adjusting for feature length. Thus, we cannot sum counts over a set of 
#' features to get the expression of that set (for example, we cannot sum counts
#'  over transcripts to get the expression for a gene). See the following link 
#'  for a discussion of RNA-seq expression units by Harold Pimentel:
#' \url{https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/}.
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' object <- summariseExprsAcrossFeatures(object)
#' }
#' 
summariseExprsAcrossFeatures <- function(object, use_as_exprs="tpm", 
                                         summarise_by="feature_id") {
    if( !is(object, "SCESet") )
        stop("Object must be an SCESet")
    if( !(summarise_by %in% colnames(fData(object))) )
        stop("The summarise_by argument is not a column of fData(object).")
    ## Define an expression matrix depending on which values we're using
    use_as_exprs <- match.arg(use_as_exprs, c("exprs", "tpm", "fpkm", "counts"))
    exprs_mat <- switch(use_as_exprs,
                        exprs=exprs(object),
                        tpm=tpm(object),
                        fpkm=fpkm(object),
                        counts=counts(object))
    ## Use reshape2 to make a long version of the expression matrix
    tmp_exprs <- data.frame(feature=fData(object)[[summarise_by]], exprs_mat)
    tmp_exprs_long <- reshape2::melt(tmp_exprs)
    exprs_new <- reshape2::acast(tmp_exprs_long, feature ~ variable, sum)
    cat("Collapsing counts to", nrow(exprs_new), "features.")
    ## Create a new SCESet object
    pd <- new("AnnotatedDataFrame", pData(object))
    fd <- new("AnnotatedDataFrame", data.frame(exprs_collapsed_to=rownames(exprs_new)))
    rownames(fd) <- rownames(exprs_new)
    sce_out <- switch(use_as_exprs,
                      exprs=newSCESet(exprsData=exprs_new, phenoData=pd, 
                                      featureData=fd),
                      tpm=newSCESet(tpmData=exprs_new, phenoData=pd, 
                                    featureData=fd),
                      fpkm=newSCESet(fpkmData=exprs_new, phenoData=pd, 
                                     featureData=fd),
                      counts=newSCESet(countData=exprs_new, phenoData=pd, 
                                       featureData=fd))
    ## Use gene symbols for rownames
    sce_out
}

# sce_kall_mmuse_gene <- summariseExprsAcrossFeatures(sce_kall_mmus)

# # Compare these counts to those produced by WTCHG Core
# load("~/021_Cell_Cycle/cache/scs_mouse_unfiltered_vers2.RData")
# scs_mouse
# head(featureNames(scs_mouse))
# ## Keep overlapping genes for kallisto data
# sum(colnames(counts_gene) %in% scs_mouse$Readgroup)
# keep_col <- colnames(counts_gene) %in% scs_mouse$Readgroup
# keep_row <- rownames(counts_gene) %in% featureNames(scs_mouse)
# counts_kallisto <- counts_gene[keep_row, keep_col]
# dim(counts_kallisto)
# ## Keep overlapping genes for Core data
# keep_row <- featureNames(scs_mouse) %in% rownames(counts_kallisto)
# keep_col <- scs_mouse$Readgroup %in% colnames(counts_kallisto)
# scs_mouse_compare <- scs_mouse[keep_row, keep_col]
# dim(scs_mouse_compare)
# ## Get genes and cells ordered the same way
# mcol <- match(scs_mouse_compare$Readgroup, colnames(counts_kallisto))
# identical(colnames(counts_kallisto)[mcol], scs_mouse_compare$Readgroup)
# mrow <- match(featureNames(scs_mouse_compare), rownames(counts_kallisto))
# identical(rownames(counts_kallisto)[mrow], featureNames(scs_mouse_compare))
# counts_kallisto <- counts_kallisto[mrow, mcol]
# ## Plot log-counts to compare quantification
# par(mfcol=c(2,2))
# for(i in 1:4) {
#     rsamp <- sample(nrow(counts_kallisto), 1000)
#     csamp <- sample(ncol(counts_kallisto), 30)
#     plot(x=(counts(scs_mouse_compare)[rsamp, csamp] + 0.5), 
#          y=(counts_kallisto[rsamp, csamp] + 0.5),
#          pch=21, bg=scales::alpha("gray30", 0.1), col="gray50", cex=0.7,
#          panel.first=grid(), log="xy", asp=1, xlab="Counts (+0.5) from WTCHG Core",
#          ylab="Counts (+0.5) from kallisto")
#     title("Sample from 1000 genes, 30 cells")
#     abline(0,1)
# }
# 
# 
