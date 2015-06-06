### Methods for the SCESet class

################################################################################
### constructor function for SCESet class

#' Create a new SCESet object.
#' 
#' Scater requires that all data be housed in SCESet objects. SCESet extends 
#' Bioconductor's ExpressionSet class, and the same basic interface is 
#' supported. newSCESet() expects a matrix of expression values as its first 
#' argument, with rows as features (usually genes) and columns as cells. 
#' Per-feature and per-cell metadata can be supplied with the featureData and 
#' phenoData arguments, respectively. Use of these optional arguments is 
#' strongly encouraged. The SCESet also includes a slot 'counts' to store an 
#' object containing raw count data. 
#' 
#' @param exprsData expression data matrix for an experiment
#' @param countData data matrix containing raw count expression values
#' @param tpmData matrix of class \code{"numeric"} containing 
#' transcripts-per-million (TPM) expression values
#' @param fpkmData matrix of class \code{"numeric"} containing fragments per 
#' kilobase of exon per million reads mapped (FPKM) expression values
#' @param phenoData data frame containing attributes of individual cells
#' @param featureData data frame containing attributes of features (e.g. genes)
#' @param experimentData MIAME class object containing metadata data and details
#' about the experiment and dataset.
#' @param lowerDetectionLimit the minimum expression level that constitutes true
#'  expression (defaults to zero and uses count data to determine if an 
#'  observation is expressed or not)
#' @param logged logical, if a value is supplied for the exprsData argument, are
#'  the expression values already on the log2 scale, or not?
#' @param is_exprsData matrix of class \code{"logical"}, indicating whether
#'    or not each observation is above the \code{lowerDetectionLimit}.
#' @return a new SCESet object
#'
#' @details
#' SCESet objects store a matrix of expression values. These values are 
#' typically transcripts-per-million (tpm), counts-per-million (cpm), fragments 
#' per kilobase per million mapped (FPKM) or some other output from a program 
#' that calculates expression values from RNA-Seq reads. We recommend that 
#' expression values on the log2 scale are used for the 'exprs' slot in the 
#' SCESet. For example, you may wish to store raw tpm values in the 'tpm' slot 
#' and \code{log2(tpm + 1)} values in the 'exprs' slot. However, expression 
#' values could also be values from a single cell qPCR run or some other type of
#'  assay. The newSCESet function can also accept raw count values. In this case
#'  see \code{\link{calculateTPM}} and \code{\link{calculateFPKM}} for computing
#'  TPM and FPKM expression values, respectively, from counts. The function 
#'  \code{\link[edgeR]{cpm}} from the package edgeR to can be used to compute 
#'  log2(counts-per-million), if desired.
#'  
#'  An \code{SCESet} object has to have the \code{'exprs'} slot defined, so if
#'  the \code{exprsData} argument is \code{NULL}, then this function will define
#'  \code{'exprs'} with the following order of precedence: log2(TPM + 1), if 
#'  \code{tpmData} is defined; log2(FPKM + 1) if \code{fpkmData} is defined; 
#'  otherwise log2(counts-per-million + 1) are used. Note that for most analyses
#'  counts-per-million are not recommended, and if possible transcripts-per-million
#'  should be used.
#'  
#'  In many downstream functions you will likely find it most convenient if the 
#'  \code{'exprs'} values are on the log2-scale, so this is recommended.
#'  
#'  
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset
newSCESet <- function( exprsData=NULL, 
                       countData=NULL,
                       tpmData=NULL,
                       fpkmData=NULL,
                       phenoData=NULL, 
                       featureData=NULL,
                       experimentData=NULL,
                       is_exprsData=NULL,
                       lowerDetectionLimit=0,
                       logged=FALSE)
{
    ## Check that we have some expression data
    if( is.null(exprsData) & is.null(countData) & is.null(tpmData) & is.null(fpkmData))
        stop("Require at least one of exprsData, tpmData, fpkmData or countData arguments.")
    ## Check dimensions of data matrices
    
    ## Check counts are a matrix; renames is_exprsData if not null
    if( !is.null(countData) )     
        countData <- as.matrix(countData)
    if( !is.null(is_exprsData) )
        isexprs <- is_exprsData
    
    ## If no exprsData provided define is_exprs from tpmData, fpkmData or countData
    if( is.null(exprsData) ) {
        ## Define exprs data if null
        if( !is.null(tpmData) ) {
            exprsData <- log2(tpmData + 1)
            logged <- TRUE
            message("exprs(object) (i.e. exprsData) is not defined. 
Using log2(transcripts-per-million + 1) for exprs slot. See also ?calculateTPM.")
        } else {
            if( !is.null(fpkmData) ) {
                exprsData <- log2(fpkmData + 1)
                logged <- TRUE
                message("exprs(object) (i.e. exprsData) is not defined. 
Using log2(FPKM + 1) for exprs slot. See also ?calculateFPKM and ?calculateTPM.")
            } else {
                exprsData <- edgeR::cpm.default(countData, prior.count = 1, log = TRUE)
                logged <- TRUE
                message("Generating log2(counts-per-million) from counts to use as
                expression data, with prior.count = 1. See edgeR::cpm().
                        Note that counts-per-million are not recommended for most analyses. 
                        Consider using transcripts-per-million instead. See ?calculateTPM.")
            }
        }    
        ## Define isexprs if null
        if( is.null(is_exprsData) ) {
            if( !is.null(tpmData) ) {
                isexprs <- tpmData > lowerDetectionLimit
                rownames(isexprs) <- rownames(tpmData)
                colnames(isexprs) <- colnames(tpmData)
                message(paste0("Defining 'is_exprs' using TPM data and a lower TPM threshold of ", lowerDetectionLimit))
            } else { 
                if( !is.null(fpkmData) ) {
                    isexprs <- fpkmData > lowerDetectionLimit
                    rownames(isexprs) <- rownames(fpkmData)
                    colnames(isexprs) <- colnames(fpkmData)
                    message(paste0("Defining 'is_exprs' using FPKM data and a lower FPKM threshold of ", lowerDetectionLimit))    
                } else {
                    isexprs <- countData > lowerDetectionLimit
                    rownames(isexprs) <- rownames(countData)
                    colnames(isexprs) <- colnames(countData)
                    message(paste0("Defining 'is_exprs' using count data and a lower count threshold of ", lowerDetectionLimit))
                }
            }
        }        
    } else {
        exprsData <- as.matrix(exprsData)
        if( is.null(is_exprsData) ) {
            isexprs <- exprsData > lowerDetectionLimit
            rownames(isexprs) <- rownames(exprsData)
            colnames(isexprs) <- colnames(exprsData)
            message(paste0("Defining 'is_exprs' using exprsData and a lower exprs threshold of ", lowerDetectionLimit))
        }
    }
    
    ## Generate valid phenoData and featureData if not provided
    if( is.null(phenoData) )
        phenoData <- annotatedDataFrameFrom(exprsData, byrow=FALSE)
    if( is.null(featureData) ) 
        featureData <- annotatedDataFrameFrom(exprsData, byrow=TRUE)
   
    ## Check experimentData
    expData_null <- new("MIAME",
                        name="<your name here>",
                        lab="<your lab here>",
                        contact="<email address>",
                        title="<title for this dataset>",
                        abstract="An SCESet",
                        url="<your website here>",
                        other=list(
                            notes="This dataset created from ...",
                            coauthors=c("")
                        ))
    if( !is.null( experimentData ) ) {
        if( is(experimentData, "MIAME") )
            expData <- experimentData
        else {
            expData <- expData_null
            warning("experimentData supplied is not an 'MIAME' object. Thus, experimentData is being set to an empty MIAME object.\n Please supply a valid 'MIAME' class object containing experiment data to experimentData(object).")
        }
    } else {
        expData <- expData_null
    }   
    
    ## Generate new SCESet object
    assaydata <- assayDataNew("environment", exprs=exprsData, is_exprs=isexprs) 
    sceset <- new( "SCESet",
                   assayData=assaydata,
                   phenoData=phenoData, 
                   featureData=featureData, 
                   experimentData=expData,
                   lowerDetectionLimit=lowerDetectionLimit,
                   logged=logged)
    
    ## Add non-null slots to assayData for SCESet object, omitting null slots
    if( !is.null(tpmData) )
        tpm(sceset) <- tpmData
    if( !is.null(fpkmData) )
        fpkm(sceset) <- fpkmData
    if( !is.null(countData) )
        counts(sceset) <- countData
    
    ## Check validity of object    
    validObject(sceset)
    sceset
}


################################################################################
### Define validity check for SCESet class object

setValidity("SCESet", function(object) {
    ## Check that the dimensions and names of the bootstraps slot are sensible
    if( (length(object@bootstraps) != 0) && (nrow(object@bootstraps) != nrow(object)) )
        return(FALSE)
    if( (length(object@bootstraps) != 0) && (ncol(object@bootstraps) != ncol(object)) )
        return(FALSE)
    if(  (length(object@bootstraps) != 0) && 
         !identical(rownames(object@bootstraps), featureNames(object)) )
        return(FALSE)
    if(  (length(object@bootstraps) != 0) && 
         !identical(colnames(object@bootstraps), sampleNames(object)) )
        return(FALSE)
    ## Check that the dimensions of the reducedDimension slot are sensible
    if( (nrow(object@reducedDimension) != 0) && 
        (nrow(object@reducedDimension) != ncol(object)) )
        return(FALSE)
    if( (nrow(object@reducedDimension) != 0) && 
        !identical(rownames(object@reducedDimension), sampleNames(object)) )
        return(FALSE)
    ## Check that the dimensions of the cellPairwiseDistances slot are sensible
    if( (nrow(object@cellPairwiseDistances) != ncol(object@cellPairwiseDistances)) ||
        (nrow(object@cellPairwiseDistances) != 0 && 
         nrow(object@cellPairwiseDistances) != ncol(object)) )
        return(FALSE)
    if( (nrow(object@cellPairwiseDistances) != 0) && 
        (!identical(rownames(object@cellPairwiseDistances), 
                    colnames(object@cellPairwiseDistances)) ||
        !identical(rownames(object@cellPairwiseDistances), sampleNames(object))) )
        return(FALSE)
    ## Check that the dimensions of the featurePairwiseDistances slot are sensible
    if( (nrow(object@featurePairwiseDistances) != 
         ncol(object@featurePairwiseDistances)) ||
        (nrow(object@featurePairwiseDistances) != 0 && 
         nrow(object@featurePairwiseDistances) != nrow(object)) )
        return(FALSE)
    if( (nrow(object@featurePairwiseDistances) != 0) && 
        (!identical(rownames(object@featurePairwiseDistances), 
                    colnames(object@featurePairwiseDistances)) ||
         !identical(rownames(object@featurePairwiseDistances), featureNames(object))) )
        return(FALSE)
    ## Check that we have sensible values for the counts
    if( is.null(counts(object)) )
        return(TRUE)
    else {
        if( any(counts(object) < 0) ) {
            warning( "The count data contain negative values." )
            return(FALSE)         
        } else
            return(TRUE)
    }
})


################################################################################
### subsetting an SCESet object

#' Subsetting SCESet Objects
#'
#' Subset method for SCESet objects, which subsets both the expression data, 
#' phenotype data, feature data and other slots in the object.
#'
#' @rdname SCESet-subset
#' @name SCESet-subset
NULL

#' @inheritParams base::Extract
#' @aliases [,SCESet,ANY-method [,SCESet,ANY,ANY-method [,SCESet,ANY,ANY,ANY-method 
#' @rdname SCESet-subset
#' @export
#' @seealso \code{\link[base]{Extract}}
#' 
setMethod('[', 'SCESet', function (x, i, j, ..., drop=FALSE) {
    if( !missing(i) && missing(j) ){
        ## Subsetting features only
        x <- selectMethod('[', 'ExpressionSet')(x, i, , drop=drop)
        if( nrow(x@featurePairwiseDistances) != 0 )
            x@featurePairwiseDistances <- x@featurePairwiseDistances[i, i, drop=drop]
        if( nrow(x@bootstraps) != 0 )
            x@bootstraps <- x@bootstraps[i, , ..., drop=drop]
    } else if( missing(i) && !missing(j) ){
        ## Subsetting cells only
        x <- selectMethod('[', 'ExpressionSet')(x, , j, drop=drop)
        if( nrow(x@cellPairwiseDistances) != 0 )
            x@cellPairwiseDistances <- x@cellPairwiseDistances[j, j, drop=drop]
        if( nrow(x@reducedDimension) != 0 )
            x@reducedDimension <- x@reducedDimension[j, , drop=drop]
        if( ncol(x@bootstraps) != 0 )
            x@bootstraps <- x@bootstraps[, j, ..., drop=drop]
    } else if( !missing(i) && !missing(j) ){
        ## Subsetting features (i) and cells (j)
        x <- selectMethod('[', 'ExpressionSet')(x, i, j, drop=drop)
        if( nrow(x@featurePairwiseDistances) != 0 )
            x@featurePairwiseDistances <- x@featurePairwiseDistances[i, i, drop=drop]
        if( nrow(x@cellPairwiseDistances) != 0 )
            x@cellPairwiseDistances <- x@cellPairwiseDistances[j, j, drop=drop]
        if( nrow(x@reducedDimension) != 0 )
            x@reducedDimension <- x@reducedDimension[j, , drop=drop]
        if( nrow(x@bootstraps) != 0 )
            x@bootstraps <- x@bootstraps[i, j, ..., drop=drop]
    } else{ 
        ## All missing: possibly not missing k for subsetting bootstrap samples
        if( nrow(x@bootstraps) != 0 && ncol(x@bootstraps) != 0 )
            x@bootstraps <- x@bootstraps[, , ..., drop=drop]
    }
    ## Check validity of object    
    validObject(x)
    x
})


################################################################################
## cellNames

#' Get cell names from an SCESet object
#' 
#' @param object An SCESet object.
#' 
#' @return A vector of cell names.
#' 
#' @details Simply a wrapper to \code{\link[Biobase]{sampleNames}}.
#' 
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' cellNames(example_sceset)
#' 
cellNames <- function(object) {
    sampleNames(object)
}


################################################################################
### Replacer methods for slots in an SCESet object

#' Replaces featureData in an SCESet object
#'
#' SCESet objects contain feature information (inherited from the ExpressionSet
#' class). This function conveniently replaces the feature data with the 
#' value supplied, which must be an AnnotatedDataFrame.
#' @param x An SCESet object.
#' @param value an AnnotatedDataFrame with updated featureData to replace 
#' existing 
#' @return A matrix of expression count data, where rows correspond to features
#' (e.g. genes) and columns correspond to cells.
#' 
#' @export
#' @rdname fData
#' @aliases fData fData,SCESet-method fData<-,SCESet,AnnotatedDataFrame-method fData<-,SCESet,data.frame-method 
#' 
#' @examples
#' \dontrun{
#' 
#' }
setReplaceMethod("fData", signature(x = "SCESet", value = "AnnotatedDataFrame"), 
                 function(x, value) {
                     x@featureData <- value
                     x
                 } )

#' @name fData
#' @rdname fData
#' @exportMethod "fData<-"
setReplaceMethod("fData", signature(x = "SCESet", value = "data.frame"), 
                 function(x, value) {
                     x@featureData <- new("AnnotatedDataFrame", value)
                     x
                 } )



#' Replaces phenoData in an SCESet object
#'
#' SCESet objects contain phenotype information (inherited from the 
#' ExpressionSet class). This function conveniently replaces the phenotype data 
#' with the value supplied, which must be an AnnotatedDataFrame.
#' @param x An SCESet object.
#' @param value an AnnotatedDataFrame with updated phenoData to replace 
#' existing 
#' @return A matrix of expression count data, where rows correspond to features
#' (e.g. genes) and columns correspond to cells.
#' 
#' @exportMethod "pData<-"
#' @rdname pData
#' @aliases pData pData,SCESet-method pData<-,SCESet,AnnotatedDataFrame-method pData<-,SCESet,data.frame-method
#' 
#' @examples
#' \dontrun{
#' 
#' }
setReplaceMethod("pData", signature(x = "SCESet", value = "AnnotatedDataFrame"), 
                 function(x, value) {
                     x@phenoData <- value
                     x
                 } )


#' @name pData
#' @rdname pData
#' @exportMethod "pData<-"
setReplaceMethod("pData", signature(x = "SCESet", value = "data.frame"), 
                 function(x, value) {
                     x@phenoData <- new("AnnotatedDataFrame", value)
                     x
                 } )


################################################################################
### counts

#' Accessors for the 'counts' element of an SCESet object.
#'
#' The counts element holds the count data as a matrix of non-negative integer
#' count values, one row for each feature (gene, exon, region, etc), and one
#' column for each cell. It is an element of the assayData slot of the SCESet
#' object.
#'
#' @usage
#' \S4method{counts}{SCESet}(object)
#'
#' \S4method{counts}{SCESet,matrix}(object)<-value
#'
#' @docType methods
#' @name counts
#' @rdname counts
#' @aliases counts counts,SCESet-method counts<-,SCESet,matrix-method
#'
#' @param object a \code{SCESet} object.
#' @param value an integer matrix
#' @author Davis McCarthy
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData = sc_example_counts)
#' counts(example_sceset)
#'
counts.SCESet <- function(object) {
    object@assayData$counts
}

#' @rdname counts
#' @export
setMethod("counts", signature(object="SCESet"), counts.SCESet)

#' @name counts
#' @rdname counts
#' @exportMethod "counts<-"
setReplaceMethod("counts", signature(object="SCESet", value="matrix"),
                 function( object, value ) {
                     object@assayData$counts <- value
                     validObject(object)
                     object
                 })

################################################################################
### is_exprs

#' Accessors for the 'is_exprs' element of an SCESet object.
#'
#' The is_exprs element holds a logical matrix indicating whether or not each 
#' observation is above the defined lowerDetectionLimit in the SCESet object. It
#' has the same dimensions as the 'exprs' and 'counts' elements, which hold the 
#' transformed expression data and count data, respectively.
#' 
#' @usage
#' \S4method{is_exprs}{SCESet}(object)
#'
#' \S4method{is_exprs}{SCESet,matrix}(object)<-value
#'
#' @docType methods
#' @name is_exprs
#' @rdname is_exprs
#' @aliases is_exprs is_exprs,SCESet-method is_exprs<-,SCESet,matrix-method
#'
#' @param object a \code{SCESet} object.
#' @param value an integer matrix
#' @author Davis McCarthy
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData = sc_example_counts)
#' is_exprs(example_sceset)
#'
is_exprs.SCESet <- function(object) {
    object@assayData$is_exprs
}

#' @rdname is_exprs
#' @export
setMethod("is_exprs", signature(object = "SCESet"), is_exprs.SCESet)

#' @name is_exprs<-
#' @rdname is_exprs
#' @exportMethod "is_exprs<-"
setReplaceMethod("is_exprs", signature(object="SCESet", value="matrix"),
                 function(object, value) {
                     object@assayData$is_exprs <- value
                     validObject(object)
                     object
                 })

################################################################################
### norm_exprs

#' Accessors for the 'norm_exprs' (normalised expression) element of an SCESet object.
#'
#' The \code{norm_exprs} element of the arrayData slot in an SCESet object holds
#' a matrix containing normalised expression values. It has the same dimensions 
#' as the 'exprs' and 'counts' elements, which hold the transformed expression 
#' data and count data, respectively.
#' 
#' @usage
#' \S4method{norm_exprs}{SCESet}(object)
#'
#' \S4method{norm_exprs}{SCESet,matrix}(object)<-value
#'
#' @docType methods
#' @name norm_exprs
#' @rdname norm_exprs
#' @aliases norm_exprs norm_exprs,SCESet-method norm_exprs<-,SCESet,matrix-method
#'
#' @param object a \code{SCESet} object.
#' @param value an integer matrix
#' 
#' @details The default for normalised expression values is mean-centred and 
#' variance-standardised expression data from the \code{exprs} slot of the 
#' \code{SCESet} object. The function \code{normaliseExprs} (or 
#' \code{normalizeExprs}) provides more options and functionality for 
#' normalising expression data.
#' 
#' @author Davis McCarthy
#' @export
#' @aliases norm_exprs norm_exprs,SCESet-method, norm_exprs<-,SCESet,matrix-method 
#' 
#' 
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData = sc_example_counts)
#' norm_exprs(example_sceset)
#'
norm_exprs.SCESet <- function(object) {
    object@assayData$norm_exprs
}

#' @rdname norm_exprs
#' @export
setMethod("norm_exprs", signature(object="SCESet"), norm_exprs.SCESet)

#' @name norm_exprs<-
#' @rdname norm_exprs
#' @exportMethod "norm_exprs<-"
setReplaceMethod("norm_exprs", signature(object="SCESet", value="matrix"),
                 function(object, value) {
                     object@assayData$norm_exprs <- value
                     validObject(object)
                     object
                 })


################################################################################
### stand_exprs

#' Accessors for the 'stand_exprs' (standardised expression) element of an SCESet object.
#'
#' The \code{stand_exprs} element of the arrayData slot in an SCESet object holds
#' a matrix containing standardised (mean-centred, variance standardised, by 
#' feature) expression values. It has the same dimensions as the 'exprs' and 
#' 'counts' elements, which hold the transformed expression data and count data,
#'  respectively.
#' 
#' @usage
#' \S4method{stand_exprs}{SCESet}(object)
#'
#' \S4method{stand_exprs}{SCESet,matrix}(object)<-value
#'
#' @docType methods
#' @name stand_exprs
#' @rdname stand_exprs
#' @aliases stand_exprs stand_exprs,SCESet-method stand_exprs<-,SCESet,matrix-method
#'
#' @param object a \code{SCESet} object.
#' @param value an integer matrix
#' 
#' @details The default for normalised expression values is mean-centred and 
#' variance-standardised expression data from the \code{exprs} slot of the 
#' \code{SCESet} object. The function \code{normaliseExprs} (or 
#' \code{normalizeExprs}) provides more options and functionality for 
#' normalising expression data.
#' 
#' @author Davis McCarthy
#' @export
#' @aliases stand_exprs stand_exprs,SCESet-method, stand_exprs<-,SCESet,matrix-method 
#' 
#' 
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData = sc_example_counts)
#' stand_exprs(example_sceset)
#'
stand_exprs.SCESet <- function(object) {
    object@assayData$stand_exprs
}

#' @rdname stand_exprs
#' @export
setMethod("stand_exprs", signature(object="SCESet"), stand_exprs.SCESet)

#' @name stand_exprs<-
#' @rdname stand_exprs
#' @exportMethod "stand_exprs<-"
setReplaceMethod("stand_exprs", signature(object="SCESet", value="matrix"),
                 function(object, value) {
                     object@assayData$stand_exprs <- value
                     validObject(object)
                     object
                 })



################################################################################
### tpm

#' Accessors for the 'tpm' (transcripts per million) element of an SCESet object.
#'
#' The \code{tpm} element of the arrayData slot in an SCESet object holds
#' a matrix containing transcripts-per-million values. It has the same dimensions 
#' as the 'exprs' and 'counts' elements, which hold the transformed expression 
#' data and count data, respectively.
#' 
#' @usage
#' \S4method{tpm}{SCESet}(object)
#'
#' \S4method{tpm}{SCESet,matrix}(object)<-value
#'
#' @docType methods
#' @name tpm
#' @rdname tpm
#' @aliases tpm tpm,SCESet-method tpm<-,SCESet,matrix-method
#'
#' @param object a \code{SCESet} object.
#' @param value a matrix of class \code{"numeric"}
#' 
#' @author Davis McCarthy
#' @export
#' @aliases tpm tpm,SCESet-method tpm<-,SCESet,matrix-method
#' 
#' @examples
#' \dontrun{
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData = sc_example_counts)
#' tpm(example_sceset)
#' }
tpm.SCESet <- function(object) {
    object@assayData$tpm
}

#' @name tpm
#' @rdname tpm
#' @export
#' @aliases tpm,SCESet-method
setMethod("tpm", signature(object="SCESet"), tpm.SCESet)

#' @name tpm<-
#' @rdname tpm
#' @exportMethod "tpm<-"
#' @aliases tpm<-,SCESet,matrix-method
setReplaceMethod("tpm", signature(object="SCESet", value="matrix"),
                 function(object, value) {
                     object@assayData$tpm <- value
                     validObject(object)
                     object
                 })


################################################################################
### cpm

#' Accessors for the 'cpm' (counts per million) element of an SCESet object.
#'
#' The \code{cpm} element of the arrayData slot in an SCESet object holds
#' a matrix containing counts-per-million values. It has the same dimensions 
#' as the 'exprs' and 'counts' elements, which hold the transformed expression 
#' data and count data, respectively.
#' 
#' @usage
#' \S4method{cpm}{SCESet}(object)
#'
#' \S4method{cpm}{SCESet,matrix}(object)<-value
#'
#' @docType methods
#' @name cpm
#' @rdname cpm
#' @aliases cpm cpm,SCESet-method cpm<-,SCESet,matrix-method
#'
#' @param object a \code{SCESet} object.
#' @param value a matrix of class \code{"numeric"}
#' 
#' @author Davis McCarthy
#' @export
#' @aliases cpm cpm,SCESet-method cpm<-,SCESet,matrix-method
#' 
#' @examples
#' \dontrun{
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData=sc_example_counts)
#' cpm(example_sceset)
#' }
cpm.SCESet <- function(object) {
    object@assayData$cpm
}

#' @name cpm
#' @rdname cpm
#' @export
#' @aliases cpm,SCESet-method
setMethod("cpm", signature(object="SCESet"), cpm.SCESet)

#' @name cpm<-
#' @rdname cpm
#' @exportMethod "cpm<-"
#' @aliases cpm<-,SCESet,matrix-method
setReplaceMethod("cpm", signature(object="SCESet", value="matrix"),
                 function(object, value) {
                     object@assayData$cpm <- value
                     validObject(object)
                     object
                 })


################################################################################
### fpkm

#' Accessors for the 'fpkm' (fragments per kilobase of exon per million reads mapped) element of an SCESet object.
#'
#' The \code{fpkm} element of the arrayData slot in an SCESet object holds
#' a matrix containing fragments per kilobase of exon per million reads mapped 
#' (FPKM) values. It has the same dimensions as the 'exprs' and 'counts' 
#' elements, which hold the transformed expression data and count data, 
#' respectively.
#' 
#' @usage
#' \S4method{fpkm}{SCESet}(object)
#'
#' \S4method{fpkm}{SCESet,matrix}(object)<-value
#'
#' @docType methods
#' @name fpkm
#' @rdname fpkm
#' @aliases fpkm fpkm,SCESet-method fpkm<-,SCESet,matrix-method
#'
#' @param object a \code{SCESet} object.
#' @param value a matrix of class \code{"numeric"}
#' 
#' @author Davis McCarthy
#' @export
#' 
#' @examples
#' \dontrun{
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData = sc_example_counts)
#' fpkm(example_sceset)
#' }
fpkm.SCESet <- function(object) {
    object@assayData$fpkm
}

#' @name fpkm
#' @rdname fpkm
#' @export
#' @aliases fpkm,SCESet-method
setMethod("fpkm", signature(object="SCESet"), fpkm.SCESet)

#' @name fpkm<-
#' @rdname fpkm
#' @exportMethod "fpkm<-"
#' @aliases fpkm<-,SCESet,matrix-method
setReplaceMethod("fpkm", signature(object="SCESet", value="matrix"),
                 function(object, value) {
                     object@assayData$fpkm <- value
                     validObject(object)
                     object
                 })




################################################################################
### bootstraps

#' Accessor and replacement for bootstrap results in an SCESet object
#'
#' SCESet objects can contain an of bootstrap expression values (for example, as 
#' generated by the kallisto software for quantifying feature abundance). These
#'  functions conveniently access and replace the 'bootstrap' slot with the value 
#'  supplied, which must be an matrix of the correct size, namely the same 
#'  number of rows and columns as the \code{SCEset} object as a whole. 
#' 
#' @docType methods
#' @name bootstraps
#' @rdname bootstraps
#' @aliases bootstraps bootstraps,SCESet-method bootstraps<-,SCESet,array-method
#'
#' @param object a \code{SCESet} object.
#' @param value an array of class \code{"numeric"} containing bootstrap 
#' expression values
#' @author Davis McCarthy
#' 
#' @return If accessing bootstraps slot of an \code{SCESet}, then an array with 
#' the bootstrap values, otherwise an \code{SCESet} object containing new 
#' bootstrap values.
#' 
#' @export
#' @aliases bootstraps bootstraps,SCESet-method bootstraps<-,SCE-Set,array-method
#' 
#' @examples
#' \dontrun{
#' ## If 'object' is an SCESet:
#' bootstraps(object)
#' }
#' 
bootstraps.SCESet <- function(object) {
    object@bootstraps
}

#' @rdname bootstraps
#' @aliases bootstraps
#' @export
setMethod("bootstraps", signature(object="SCESet"), bootstraps.SCESet)


#' @name bootstraps<-
#' @aliases bootstraps
#' @rdname bootstraps
#' @export "bootstraps<-"
setReplaceMethod("bootstraps", signature(object="SCESet", value="array"), 
                 function(object, value) {
                     if( (nrow(value) == nrow(object)) & (ncol(value) == ncol(object)) ) {
                         object@bootstraps <- value
                         return(object)
                     } else
                         stop("Array supplied is of incorrect size.")
                 } )


################################################################################
### reducedDimension

#' Reduced dimension representation for cells in an SCESet object
#'
#' SCESet objects can contain a matrix of reduced dimension coordinates for 
#' cells. These functions conveniently access and replace the reduced dimension 
#' coordinates with the value supplied, which must be a matrix of the correct 
#' size. The function \code{redDim} is simply shorthand for 
#' \code{reducedDimension}.
#' 
#' @docType methods
#' @name reducedDimension
#' @rdname reducedDimension
#' @aliases reducedDimension reducedDimension,SCESet-method reducedDimension<-,SCESet,matrix-method redDim,SCESet-method redDim<-,SCESet,matrix-method
#'
#' @param object a \code{SCESet} object.
#' @param value a matrix of class \code{"numeric"} containing reduced dimension
#' coordinates for cells.
#' @author Davis McCarthy
#' 
#' @return If accessing the \code{reducedDimension} slot, then the matrix of 
#' reduced dimension coordinates. If replacing the \code{reducedDimension} slot
#' then the new matrix is added to the \code{SCESet} object.
#' 
#' @export
#' @examples
#' \dontrun{
#' 
#' }
#' 
reducedDimension.SCESet <- function(object) {
    object@reducedDimension
}

#' @rdname reducedDimension
#' @aliases reducedDimension
#' @export
setMethod("reducedDimension", signature(object="SCESet"), 
          reducedDimension.SCESet)

#' @rdname reducedDimension
#' @aliases reducedDimension
#' @export
redDim.SCESet <- function(object) {
    object@reducedDimension
}

#' @rdname reducedDimension
#' @aliases reducedDimension
#' @export
setMethod("redDim", signature(object="SCESet"), redDim.SCESet)

#' @name reducedDimension<-
#' @aliases reducedDimension
#' @rdname reducedDimension
#' @exportMethod "reducedDimension<-"
setReplaceMethod("reducedDimension", signature(object="SCESet", value="matrix"), 
                 function(object, value) {
                     if( nrow(value) == ncol(object) ) {
                         object@reducedDimension <- value
                         return(object)
                     }
                     else
                         stop("Reduced dimension matrix supplied is of incorrect size. 
                              Rows of reduced dimension matrix should correspond to cells, i.e. columns of SCESet object.")
                 } )

#' @name redDim<-
#' @aliases reducedDimension
#' @rdname reducedDimension
#' @exportMethod "redDim<-"
setReplaceMethod("redDim", signature(object="SCESet", value="matrix"), 
                 function(object, value) {
                     if( nrow(value) == ncol(object) ) {
                         object@reducedDimension <- value
                         return(object)
                     }
                     else
                         stop("Reduced dimension matrix supplied is of incorrect size. 
                              Rows of reduced dimension matrix should correspond to cells, i.e. columns of SCESet object.")
                 } )



################################################################################
### cellPairwiseDistances

#' cellPairwiseDistances in an SCESet object
#'
#' SCESet objects can contain a matrix of pairwise distances between cells. These
#'  functions conveniently access and replace the cell pairwise distances with the value 
#'  supplied, which must be a matrix of the correct size. The function \code{cellDist}
#'  is simply shorthand for \code{cellPairwiseDistances}.
#' 
#' @docType methods
#' @name cellPairwiseDistances
#' @rdname cellPairwiseDistances
#' @aliases cellPairwiseDistances cellPairwiseDistances,SCESet-method cellPairwiseDistances<-,SCESet,matrix-method cellDist,SCESet-method cellDist<-,SCESet,matrix-method
#'
#' @param object a \code{SCESet} object.
#' @param value a matrix of class \code{"numeric"} containing cell pairwise 
#' distances
#' @author Davis McCarthy
#' 
#' @return An SCESet object containing new cell pairwise distances matrix.
#' 
#' @export
#' @examples
#' \dontrun{
#' 
#' }
#' 
cellPairwiseDistances.SCESet <- function(object) {
    object@cellPairwiseDistances
}

#' @rdname cellPairwiseDistances
#' @aliases cellPairwiseDistances
#' @export
setMethod("cellPairwiseDistances", signature(object="SCESet"), 
          cellPairwiseDistances.SCESet)

#' @rdname cellPairwiseDistances
#' @aliases cellPairwiseDistances
#' @export
cellDist.SCESet <- function(object) {
    object@cellPairwiseDistances
}

#' @rdname cellPairwiseDistances
#' @aliases cellPairwiseDistances
#' @export
setMethod("cellDist", signature(object="SCESet"), cellDist.SCESet)

#' @name cellPairwiseDistances<-
#' @aliases cellPairwiseDistances
#' @rdname cellPairwiseDistances
#' @exportMethod "cellPairwiseDistances<-"
setReplaceMethod("cellPairwiseDistances", signature(object="SCESet", value="matrix"), 
                 function(object, value) {
                     if( nrow(value) == ncol(object) ) {
                         object@cellPairwiseDistances <- value
                         return(object)
                     }
                     else
                         stop("Cell pairwise distance matrix supplied is of incorrect size.")
                 } )

#' @name cellDist<-
#' @aliases cellPairwiseDistances
#' @rdname cellPairwiseDistances
#' @exportMethod "cellDist<-"
setReplaceMethod("cellDist", signature(object="SCESet", value="matrix"), 
                 function(object, value) {
                     if( nrow(value) == ncol(object) ) {
                         object@cellPairwiseDistances <- value
                         return(object)
                     }
                     else
                         stop("Cell pairwise distance matrix supplied is of incorrect size.")
                 } )


################################################################################
### featurePairwiseDistances

#' featurePairwiseDistances in an SCESet object
#'
#' SCESet objects can contain a matrix of pairwise distances between features 
#' (e.g. genes, transcripts). These functions conveniently access and replace 
#' the gene pairwise distances with the value supplied, which must be a matrix 
#' of the correct size. The function \code{featDist} is simply shorthand for 
#' \code{featurePairwiseDistances}.
#'
#' @param object a \code{SCESet} object.
#' @param value a matrix of class \code{"numeric"} containing feature pairwise 
#' distances
#'
#' @docType methods
#' @name featurePairwiseDistances
#' @rdname featurePairwiseDistances
#' @aliases featurePairwiseDistances featurePairwiseDistances,SCESet-method featurePairwiseDistances<-,SCESet,matrix-method featDist featDist,SCESet-method featDist<-,SCESet,matrix-method
#'
#' @author Davis McCarthy
#' 
#' @return An SCESet object containing new feature pairwise distances matrix.
#' @export
#' @examples
#' \dontrun{
#' 
#' }
#' 
featurePairwiseDistances.SCESet <- function(object) {
    object@featurePairwiseDistances
}

#' @rdname featurePairwiseDistances
#' @aliases featurePairwiseDistances 
#' @export
setMethod("featurePairwiseDistances", signature(object="SCESet"), 
          featurePairwiseDistances.SCESet)

#' @aliases featurePairwiseDistances
#' @rdname featurePairwiseDistances
#' @export
featDist.SCESet <- function(object) {
    object@featurePairwiseDistances
}

#' @aliases featurePairwiseDistances
#' @rdname featurePairwiseDistances
#' @export
setMethod("featDist", signature(object="SCESet"), featDist.SCESet)

#' @name featurePairwiseDistances<-
#' @aliases featurePairwiseDistances
#' @rdname featurePairwiseDistances
#' @export "featurePairwiseDistances<-"
setReplaceMethod("featurePairwiseDistances", signature(object="SCESet", value="matrix"), 
                 function(object, value) {
                     if( nrow(value) == nrow(object) ) {
                         object@featurePairwiseDistances <- value
                         return(object)
                     }
                     else
                         stop("Feature pairwise distance matrix supplied is of incorrect size.")
                 } )

#' @name featDist<-
#' @rdname featurePairwiseDistances
#' @aliases featurePairwiseDistances
#' @export "featDist<-"
setReplaceMethod("featDist", signature(object="SCESet", value="matrix"), 
                 function(object, value) {
                     if( nrow(value) == nrow(object) ) {
                         object@featurePairwiseDistances <- value
                         return(object)
                     }
                     else
                         stop("Feature pairwise distance matrix supplied is of incorrect size.")
                 } )


################################################################################
### Convert to and from Monocle CellDataSet objects

#' Convert an \code{SCESet} to a \code{CellDataSet}
#' 
#' @param sce An \code{SCESet} object 
#' @param useExpression If TRUE (default), `exprs(sce)` is used as the `exprsData`, otherwise `counts(sce)`
#' 
#' @export
#' @importClassesFrom monocle CellDataSet
#' @rdname toCellDataSet
#' @name toCellDataSet
#' @return An object of class \code{SCESet}
toCellDataSet <- function(sce, useExpression=TRUE) {
    pkgAvail <- requireNamespace("monocle", character.only=TRUE) 
    if(pkgAvail) {
        if(!is(sce,'SCESet')) stop('sce must be of type SCESet')
        exprsData <- NULL
        if(useExpression) {
            exprsData <- exprs(sce)
        } else {
            exprsData <- counts(sce)
        }
        
        cds <- monocle::newCellDataSet(exprsData, phenoData=phenoData(sce),
                                       featureData=featureData(sce),
                                       lowerDetectionLimit=sce@lowerDetectionLimit)
        return( cds )
    } 
    else 
        stop("Require package monocle to be installed to use this function.")
}

#' Convert a \code{CellDataSet} to an \code{SCESet}
#' 
#' @param cds A \code{CellDataSet} from the \code{monocle} package
#' @param useExpression logical If TRUE (default), `exprsData` is mapped to `exprs(sce)`, otherwise `counts(sce)`
#' @param logged logical, if a value is supplied for the exprsData argument, are the expression values already on the log2 scale, or not?
#' 
#' @export
#' @importClassesFrom monocle CellDataSet
#' @rdname fromCellDataSet
#' @name fromCellDataSet
#' @return An object of class \code{SCESet}
fromCellDataSet <- function(cds, useExpression=TRUE, logged=FALSE) {
    pkgAvail <- requireNamespace("monocle", character.only=TRUE) 
    if(pkgAvail) {
        if(!is(cds,'CellDataSet')) stop('cds must be of type CellDataSet from package monocle')
        exprsData <- countData <- NULL
        if(useExpression) {
            exprsData <- exprs(cds)
        } else {
            countData <- exprs(cds)
        }
        
        sce <- newSCESet(exprsData=exprsData, phenoData=phenoData(HSMM),
                         featureData=featureData(cds), countData=countData,
                         lowerDetectionLimit=cds@lowerDetectionLimit,
                         logged=logged)
        return( sce )    
    }
    else 
        stop("Require package monocle to be installed to use this function.")
}

