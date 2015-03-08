### Methods for the SCESet class

################################################################################
### constructor function for SCESet class

#' Creates a new SCESet object.
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
#' @param cellData expression data matrix for an experiment
#' @param phenoData data frame containing attributes of individual cells
#' @param featureData data frame containing attributes of features (e.g. genes)
#' @param countData data matrix containing raw count expression values
#' @param lowerDetectionLimit the minimum expression level that constitutes true
#'  expression (defaults to zero and uses count data to determine if an 
#'  observation is expressed or not)
#'  @param logged logical, if a value is supplied for the cellData argument, are
#'  the expression values already on the log2 scale, or not?
#'  @param isExpr matrix of class \code{"logical"}, indicating whether
#'    or not each observation is above the \code{lowerDetectionLimit}.
#' @return a new SCESet object
#' @details
#' SCESet objects store a matrix of expression values. These values are 
#' typically counts-per-million (cpm), fragments per kilobase per million mapped
#' (FPKM) or some other output from a program that calculates expression values 
#' from RNA-Seq reads. However, they could also be values from a single cell 
#' qPCR run or some other type of assay. The newSCESet function can also accept 
#' raw count values, in which case it uses a function from the package edgeR to 
#' compute log2(counts-per-million) and uses these as the expression values.
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset
newSCESet <- function( cellData = NULL, 
                       phenoData = NULL, 
                       featureData = NULL,
                       countData = NULL,
                       lowerDetectionLimit = 0,
                       logged = FALSE,
                       isExpr = NULL)
{
    ## Check that we have some expression data
    if(is.null(cellData) & is.null(countData))
        stop("Require at least one of cellData or countData object")
    ## Check counts are a matrix
    if(is.null(countData)) {
        ## Have to insert NAs for counts
        countData <- matrix(NA, nrow = nrow(cellData), ncol = ncol(cellData))
        rownames(countData) <- rownames(cellData)
        colnames(countData) <- colnames(cellData)
    } else            
        countData <- as.matrix(countData)
    ## If no cellData provided, generate as cpm from count matrix
    if(is.null(cellData)) {
        cellData <- edgeR::cpm.default(countData, prior.count = 1, log = TRUE)
        logged <- TRUE
        message("Generating log2(counts-per-million) from counts to use as 
                expression data, with prior.count = 1. See edgeR::cpm().")
        if( is.null(isExpr) ) {
            isexpr <- countData > lowerDetectionLimit
            rownames(isexpr) <- rownames(countData)
            colnames(isexpr) <- colnames(countData)
            message(paste0("Defining 'isExpr' using count data and a lower count
                           threshold of ", lowerDetectionLimit))
        }        
    } else {
        cellData <- as.matrix(cellData)
        if( is.null(isExpr) ) {
            isexpr <- cellData > lowerDetectionLimit
            rownames(isexpr) <- rownames(cellData)
            colnames(isexpr) <- colnames(cellData)
            message(paste0("Defining 'isExpr' using cellData and a lower count 
                           threshold of ", lowerDetectionLimit))
        }
    }
    ## Generate valid phenoData and featureData if not provided
    if( is.null( phenoData ) )
        phenoData <- annotatedDataFrameFrom(cellData, byrow = FALSE)
    if( is.null( featureData ) ) 
        featureData <- annotatedDataFrameFrom(cellData, byrow = TRUE)
    ## Generate new SCESet object
    sceset <- new( "SCESet",
                   assayData = assayDataNew("environment", 
                                            exprs = cellData, 
                                            counts = countData,
                                            isExpr = isexpr),
                   phenoData = phenoData, 
                   featureData = featureData, 
                   lowerDetectionLimit = lowerDetectionLimit,
                   logged = logged)
    ## Check validity of object    
    validObject(sceset)
    sceset
}


################################################################################
### Define validity check for SCESet class object

setValidity("SCESet", function(object) {
    if (all( is.na( counts(object) ) ) ) {
       return(TRUE)
    }  else {
        if ( any( counts(object) < 0 ) ) {
            warning( "The count data contain negative values" )
            return(FALSE)         
        } else
            return(TRUE)
    }
} )


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
#' @export
#' @examples
#' \dontrun{
#' 
#' }
setReplaceMethod("fData", signature(x = "SCESet", value = "AnnotatedDataFrame"), 
                 function(x, value) {
                     x@featureData <- value
                     x
                 } )


#' Replaces phenoData in an SCESet object
#'
#' SCESet objects contain phenotype information (inherited from the 
#' ExpressionSet class). This function conveniently replaces the feature data 
#' with the value supplied, which must be an AnnotatedDataFrame.
#' @param x An SCESet object.
#' @param value an AnnotatedDataFrame with updated phenoData to replace 
#' existing 
#' @return A matrix of expression count data, where rows correspond to features
#' (e.g. genes) and columns correspond to cells.
#' @export
#' @examples
#' \dontrun{
#' 
#' }
setReplaceMethod("pData", signature(x = "SCESet", value = "AnnotatedDataFrame"), 
                 function(x, value) {
                     x@phenoData <- value
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
### isExpr

#' Accessors for the 'isExpr' element of an SCESet object.
#'
#' The isExpr element holds a logical matrix indicating whether or not each 
#' observation is above the defined lowerDetectionLimit in the SCESet object. It
#' has the same dimensions as the 'exprs' and 'counts' elements, which hold the 
#' transformed expression data and count data, respectively.
#' 
#' @usage
#' \S4method{isExpr}{SCESet}(object)
#'
#' \S4method{isExpr}{SCESet,matrix}(object)<-value
#'
#' @docType methods
#' @name isExpr
#' @rdname isExpr
#' @aliases isExpr isExpr,SCESet-method isExpr<-,SCESet,matrix-method
#'
#' @param object a \code{SCESet} object.
#' @param value an integer matrix
#' @author Davis McCarthy
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData = sc_example_counts)
#' isExpr(example_sceset)
#'
isExpr.SCESet <- function(object) {
    object@assayData$isExpr
}


#' @rdname isExpr
#' @export
setMethod("isExpr", signature(object = "SCESet"), isExpr.SCESet)

#' @name isExpr<-
#' @rdname isExpr
#' @export "isExpr<-"
setReplaceMethod("isExpr", signature(object = "SCESet", value = "matrix"),
                 function( object, value ) {
                     object@assayData$isExpr <- value
                     validObject(object)
                     object
                 })

