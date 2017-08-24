################################################################################
### updating an old SCESet object

#' @rdname toSingleCellExperiment
#' @export
#' @examples
#' \dontrun{
#' updateSCESet(example_sceset)
#' }
updateSCESet <- function(object) {
    toSingleCellExperiment(object)
}

#' Convert an SCESet object to a SingleCellExperiment object
#'
#' Convert an SCESet object produced with an older version of the
#' package to a SingleCellExperiment object compatible with the current version.
#'
#' @param object an \code{\link{SCESet}} object to be updated
#'
#' @return a \code{\link{SingleCellExperiment}} object
#' 
#' @name toSingleCellExperiment
#' @rdname toSingleCellExperiment
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @import SingleCellExperiment
#' @export
#' @examples
#' \dontrun{
#' toSingleCellExperiment(example_sceset)
#' }
toSingleCellExperiment <- function(object) {
    stopifnot(methods::is(object, "SCESet"))
    new.assay <- list()
    for (x in Biobase::assayDataElementNames(object)) {
        new.assay[[x]] <- Biobase::assayDataElement(object, x)
    }
    
    new.coldata <- DataFrame(row.names = Biobase::sampleNames(object))
    pdat <- Biobase::phenoData(object)
    for (x in varLabels(pdat)) { 
        new.coldata[[x]] <- pdat[[x]]
    }
    
    new.rowdata <- DataFrame(row.names = Biobase::featureNames(object)) 
    fdat <- Biobase::featureData(object)
    for (x in varLabels(fdat)) { 
        new.rowdata[[x]] <- fdat[[x]]
    }

    sce <- SingleCellExperiment(new.assay, colData = new.coldata, 
                         rowData = new.rowdata)
        
    if (nrow(object@reducedDimension) > 0 && 
        ncol(object@reducedDimension) > 0) {
        new.reddim <- SimpleList(redDim = reducedDim(object))
        reducedDim(sce) <- new.reddim
    }

    sce
}


################################################################################

#' Create a new SCESet object
#'
#' Deprecated from scater version 1.3.29; the package now uses the 
#' \code{\link{SingleCellExperiment}} class. To convert an SCESet object to 
#' SingleCellExperiment see the \code{\link{toSingleCellExperiment}} function.
#'
newSCESet <- function() {
    stop("Deprecated from scater 1.5.2. Use SingleCellExperiment() instead.")
}



################################################################################

#' Retrieve a representation of gene expression
#'
#' Deprecated from scater version 1.3.29.
#'
#' @param object An object of type \code{SCESet}
#' @return A matrix representation of expression values.
#'
getExprs <- function(object) {
    stop("Deprecated from scater 1.3.29")
}


################################################################################
### Convert to and from Monocle CellDataSet objects

#' Convert an \code{SCESet} to a \code{CellDataSet}. 
#' 
#' Deprecated from scater version 1.5.2.
#'
#' @param sce An \code{SCESet} object
#' @param exprs_values What should \code{exprs(cds)} be mapped from in the
#' \code{SCESet}? Should be one of "exprs", "tpm", "fpkm", "counts"
#'
#' @export
#' @rdname toCellDataSet
#' @name toCellDataSet
#' @return An object of class \code{CellDataSet}
#' @examples
#' \dontrun{"Deprecated"}
toCellDataSet <- function(sce, exprs_values = "exprs") {
    stop("Deprecated from scater 1.5.2.")
}

#' Convert a \code{CellDataSet} to an \code{SCESet}
#' 
#' Deprecated from scater version 1.5.2.
#'
#' @param cds A \code{CellDataSet} from the \code{monocle} package
#' @param exprs_values What should \code{exprs(cds)} be mapped to in the 
#' \code{SCESet}? Should be one of "exprs", "tpm", "fpkm", "counts"
#' @param logged logical, if \code{exprs_values="exprs"}, are
#'  the expression values already on the log2 scale, or not?
#' @param logExprsOffset numeric, value to add prior to log-transformation.
#'
#' @export
#' @importFrom Biobase featureData
#' @importFrom Biobase phenoData
#' @rdname fromCellDataSet
#' @name fromCellDataSet
#' @return An object of class \code{SCESet}
#' @examples
#' \dontrun{"Deprecated"}
fromCellDataSet <- function(cds, exprs_values = "tpm", logged = FALSE, 
                            logExprsOffset = 1) {
    stop("Deprecated from scater 1.5.2.")
}

