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
#' @param experimentData MIAME class object containing metadata data and details
#' about the experiment and dataset.
#' @param lowerDetectionLimit the minimum expression level that constitutes true
#'  expression (defaults to zero and uses count data to determine if an 
#'  observation is expressed or not)
#' @param logged logical, if a value is supplied for the cellData argument, are
#'  the expression values already on the log2 scale, or not?
#' @param isExprs matrix of class \code{"logical"}, indicating whether
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
                       experimentData = NULL,
                       countData = NULL,
                       lowerDetectionLimit = 0,
                       logged = FALSE,
                       isExprs = NULL)
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
        if( is.null(isExprs) ) {
            isexprs <- countData > lowerDetectionLimit
            rownames(isexprs) <- rownames(countData)
            colnames(isexprs) <- colnames(countData)
            message(paste0("Defining 'isExprs' using count data and a lower count
                           threshold of ", lowerDetectionLimit))
        }        
    } else {
        cellData <- as.matrix(cellData)
        if( is.null(isExprs) ) {
            isexprs <- cellData > lowerDetectionLimit
            rownames(isexprs) <- rownames(cellData)
            colnames(isexprs) <- colnames(cellData)
            message(paste0("Defining 'isExprs' using cellData and a lower count 
                           threshold of ", lowerDetectionLimit))
        }
    }
    ## Generate valid phenoData and featureData if not provided
    if( is.null( phenoData ) )
        phenoData <- annotatedDataFrameFrom(cellData, byrow = FALSE)
    if( is.null( featureData ) ) 
        featureData <- annotatedDataFrameFrom(cellData, byrow = TRUE)
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
    sceset <- new( "SCESet",
                   assayData = assayDataNew("environment", 
                                            exprs = cellData, 
                                            counts = countData,
                                            isExprs = isexprs),
                   phenoData = phenoData, 
                   featureData = featureData, 
                   experimentData = expData,
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
### isExprs

#' Accessors for the 'isExprs' element of an SCESet object.
#'
#' The isExprs element holds a logical matrix indicating whether or not each 
#' observation is above the defined lowerDetectionLimit in the SCESet object. It
#' has the same dimensions as the 'exprs' and 'counts' elements, which hold the 
#' transformed expression data and count data, respectively.
#' 
#' @usage
#' \S4method{isExprs}{SCESet}(object)
#'
#' \S4method{isExprs}{SCESet,matrix}(object)<-value
#'
#' @docType methods
#' @name isExprs
#' @rdname isExprs
#' @aliases isExprs isExprs,SCESet-method isExprs<-,SCESet,matrix-method
#'
#' @param object a \code{SCESet} object.
#' @param value an integer matrix
#' @author Davis McCarthy
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData = sc_example_counts)
#' isExprs(example_sceset)
#'
isExprs.SCESet <- function(object) {
    object@assayData$isExprs
}

#' @rdname isExprs
#' @export
setMethod("isExprs", signature(object = "SCESet"), isExprs.SCESet)

#' @name isExprs<-
#' @rdname isExprs
#' @exportMethod "isExprs<-"
setReplaceMethod("isExprs", signature(object="SCESet", value="matrix"),
                 function(object, value) {
                     object@assayData$isExprs <- value
                     validObject(object)
                     object
                 })

################################################################################
### norm_exprs

#' Accessors for the 'norm_exprs' (normalised expression) element of an SCESet object.
#'
#' The \code{norm_exrps} element of the arrayData slot in an SCESet object holds
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
### genePairwiseDistances

#' genePairwiseDistances in an SCESet object
#'
#' SCESet objects can contain a matrix of pairwise distances between genes. These
#'  functions conveniently access and replace the gene pairwise distances with the value 
#'  supplied, which must be a matrix of the correct size. The function \code{geneDist}
#'  is simply shorthand for \code{genePairwiseDistances}.
#'
#' @param object a \code{SCESet} object.
#' @param value a matrix of class \code{"numeric"} containing gene pairwise 
#' distances
#'
#' @docType methods
#' @name genePairwiseDistances
#' @rdname genePairwiseDistances
#' @aliases genePairwiseDistances genePairwiseDistances,SCESet-method genePairwiseDistances<-,SCESet,matrix-method geneDist geneDist,SCESet-method geneDist<-,SCESet,matrix-method
#'
#' @author Davis McCarthy
#' 
#' @return An SCESet object containing new gene pairwise distances matrix.
#' @export
#' @examples
#' \dontrun{
#' 
#' }
#' 
genePairwiseDistances.SCESet <- function(object) {
    object@genePairwiseDistances
}

#' @rdname genePairwiseDistances
#' @aliases genePairwiseDistances 
#' @export
setMethod("genePairwiseDistances", signature(object="SCESet"), 
          genePairwiseDistances.SCESet)

#' @aliases genePairwiseDistances
#' @rdname genePairwiseDistances
#' @export
geneDist.SCESet <- function(object) {
    object@genePairwiseDistances
}

#' @aliases genePairwiseDistances
#' @rdname genePairwiseDistances
#' @export
setMethod("geneDist", signature(object="SCESet"), geneDist.SCESet)

#' @name genePairwiseDistances<-
#' @aliases genePairwiseDistances
#' @rdname genePairwiseDistances
#' @export "genePairwiseDistances<-"
setReplaceMethod("genePairwiseDistances", signature(object="SCESet", value="matrix"), 
                 function(object, value) {
                     if( nrow(value) == nrow(object) ) {
                         object@genePairwiseDistances <- value
                         return(object)
                     }
                     else
                         stop("Gene pairwise distance matrix supplied is of incorrect size.")
                 } )

#' @name geneDist<-
#' @rdname genePairwiseDistances
#' @aliases genePairwiseDistances
#' @export "geneDist<-"
setReplaceMethod("geneDist", signature(object="SCESet", value="matrix"), 
                 function(object, value) {
                     if( nrow(value) == nrow(object) ) {
                         object@genePairwiseDistances <- value
                         return(object)
                     }
                     else
                         stop("Gene pairwise distance matrix supplied is of incorrect size.")
                 } )


################################################################################
### Convert to and from Monocle CellDataSet objects

#' Convert an \code{SCESet} to a \code{CellDataSet}
#' 
#' @param sce An \code{SCESet} object 
#' @param useExpression If TRUE (default), `exprs(sce)` is used as the `cellData`, otherwise `counts(sce)`
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
        cellData <- NULL
        if(useExpression) {
            cellData <- exprs(sce)
        } else {
            cellData <- counts(sce)
        }
        
        cds <- monocle::newCellDataSet(cellData, phenoData=phenoData(sce),
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
#' @param useExpression logical If TRUE (default), `cellData` is mapped to `exprs(sce)`, otherwise `counts(sce)`
#' @param logged logical, if a value is supplied for the cellData argument, are the expression values already on the log2 scale, or not?
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
        cellData <- countData <- NULL
        if(useExpression) {
            cellData <- exprs(cds)
        } else {
            countData <- exprs(cds)
        }
        
        sce <- newSCESet(cellData=cellData, phenoData=phenoData(HSMM),
                         featureData=featureData(cds), countData=countData,
                         lowerDetectionLimit=cds@lowerDetectionLimit,
                         logged=logged)
        return( sce )    
    }
    else 
        stop("Require package monocle to be installed to use this function.")
}

