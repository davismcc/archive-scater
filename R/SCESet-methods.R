################################################################################
### updating an old SCESet object

#' Update an SCESet object to the current version
#'
#' It can be necessary to update an SCESet produced with an older version of the
#' package to be compatible with the current version of the package.
#'
#' @param object an \code{\link{SCESet}} object to be updated
#'
#' @return an updated \code{\link{SCESet}} object
#'
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' updateSCESet(example_sceset)
updateSCESet <- function(object) {
    new.assay <- Assays()
    for (x in assayDataElementNames(object)) {
        new.assay[[x]] <- assayDataElement(object, x)
    }

    new.coldata <- DataFrame(row.names=sampleNames(object))
    pdat <- phenoData(object)
    for (x in varLabels(pdat)) { 
        new.coldata[[x]] <- pdat[[x]]
    }
    
    new.rowdata <- DataFrame(row.names=featureNames(object)) 
    fdat <- featureData(object)
    for (x in varLbels(fdat)) { 
        new.rowdata[[x]] <- fdat[[x]]
    }

    SingleCellExperiment(new.assay, colData=new.coldata, rolData=new.rowdata)
}

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
#' @importFrom BiocGenerics counts
#' @aliases counts counts,SCESet-method counts<-,SCESet,matrix-method
#' @return A matrix of count values.
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
#' @rdname counts
#' @export
setMethod("counts", "SingleCellExperiment", function(object) {
    assay(object, i="counts")
})

#' @name counts
#' @rdname counts
#' @importFrom BiocGenerics counts<-
#' @exportMethod "counts<-"
setReplaceMethod("counts", c("SingleCellExperiment", "ANY"), function(object, value) {
    assay(object, i="counts") <- value
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


#' @rdname norm_exprs
#' @export
setMethod("norm_exprs", "SingleCellExperiment", function(object) { 
    assay(object, i="norm_exprs")
})

#' @name norm_exprs<-
#' @rdname norm_exprs
#' @exportMethod "norm_exprs<-"
setReplaceMethod("norm_exprs", c("SingleCellExperiment", "ANY"), function(object, value) {                 
    assay(object, i="norm_exprs") <- value
    object
})

################################################################################
### stand_exprs

#' Accessors for the 'stand_exprs' (standardised expression) element of an
#' SCESet object.
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

#' @rdname stand_exprs
#' @export
setMethod("stand_exprs", "SingleCellExperiment", function(object) { 
    assay(object, i="stand_exprs")
})

#' @name stand_exprs<-
#' @rdname stand_exprs
#' @exportMethod "stand_exprs<-"
setReplaceMethod("stand_exprs", c("SingleCellExperiment", "ANY"), function(object, value) {                 
    assay(object, i="stand_exprs") <- value
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
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' bootstraps(example_sceset)
#'
bootstraps.SCESet <- function(object) {
    object@bootstraps
}

#' @rdname bootstraps
#' @aliases bootstraps
#' @export
setMethod("bootstraps", signature(object = "SCESet"), bootstraps.SCESet)


#' @name bootstraps<-
#' @aliases bootstraps
#' @rdname bootstraps
#' @export "bootstraps<-"
setReplaceMethod("bootstraps", signature(object = "SCESet", value = "array"),
                 function(object, value) {
                     if ( (nrow(value) == nrow(object)) &&
                          (ncol(value) == ncol(object)) ) {
                         object@bootstraps <- value
                         return(object)
                     } else
                         stop("Array supplied is of incorrect size.")
                 } )

################################################################################
### Convert to and from Monocle CellDataSet objects

#' Convert an \code{SCESet} to a \code{CellDataSet}
#'
#' @param sce An \code{SCESet} object
#' @param exprs_values What should \code{exprs(cds)} be mapped from in the \code{SCESet}? Should be
#' one of "exprs", "tpm", "fpkm", "counts"
#'
#' @export
#' @rdname toCellDataSet
#' @name toCellDataSet
#' @return An object of class \code{SCESet}
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' if ( requireNamespace("monocle") ) {
#'     toCellDataSet(example_sceset)
#' }
toCellDataSet <- function(sce, exprs_values = "exprs") {
    pkgAvail <- requireNamespace("monocle")
    if (pkgAvail) {
        if (!is(sce,'SCESet')) stop('sce must be of type SCESet')
        exprsData <- NULL
        exprs_values <- match.arg(exprs_values, c("exprs", "tpm", "fpkm", "counts"))
        exprsData <- switch(exprs_values,
                            exprs = exprs(sce),
                            tpm = tpm(sce),
                            fpkm = fpkm(sce),
                            counts = counts(sce))
        if ( exprs_values == "exprs" ) 
            exprsData <- 2 ^ exprsData - sce@logExprsOffset
        celldataset <- monocle::newCellDataSet(
            exprsData, phenoData = phenoData(sce),
            featureData = featureData(sce),
            lowerDetectionLimit = sce@lowerDetectionLimit)
        celldataset@reducedDimS <- t(redDim(sce))
        celldataset
    }
    else
        stop("Require package monocle to be installed to use this function.")
}

#' Convert a \code{CellDataSet} to an \code{SCESet}
#'
#' @param cds A \code{CellDataSet} from the \code{monocle} package
#' @param exprs_values What should \code{exprs(cds)} be mapped to in the \code{SCESet}? Should be
#' one of "exprs", "tpm", "fpkm", "counts"
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
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' if ( requireNamespace("monocle") ) {
#'     # cds <- toCellDataSet(example_sceset) # not run
#'     # sceset <- fromCellDataSet(cds) # not run
#' }
fromCellDataSet <- function(cds, exprs_values = "tpm", logged = FALSE, logExprsOffset = 1) {
    pkgAvail <- requireNamespace("monocle")
    if (pkgAvail) {
        if (!is(cds,'CellDataSet')) stop('cds must be of type CellDataSet from package monocle')
        exprsData <- countData <- NULL
        exprs_values <- match.arg(exprs_values,
                                  c("exprs", "tpm", "fpkm", "counts"))
        exprsData <- countData <- tpmData <- fpkmData <- NULL
        if (exprs_values == "exprs") {
            exprsData <- exprs(cds)
            if (!logged) exprsData <- log2(exprsData + logExprsOffset)
        } else if (exprs_values == "tpm") {
            tpmData <- exprs(cds)
        } else if (exprs_values == "fpkm") {
            fpkmData <- exprs(cds)
        } else {
            countData <- exprs(cds)
        }

        sce <- newSCESet(exprsData = exprsData, tpmData = tpmData,
                         fpkmData = fpkmData, countData = countData,
                         phenoData = phenoData(cds),
                         featureData = featureData(cds),
                         lowerDetectionLimit = cds@lowerDetectionLimit,
                         logExprsOffset = logExprsOffset)

        ## now try and preserve a reduced dimension representation
        ## this is really not elegant - KC
        rds <- cds@reducedDimS
        if (length(rds) == 0) {
            rds <- cds@reducedDimA
        }
        if (length(rds) == 2 * ncol(sce)) { # something is there and of the right dimension
            if (nrow(rds) == 2) {
                redDim(sce) <- t(rds)
            } else if (ncol(rds) == 2) {
                redDim(sce) <- rds
            } # else do nothing
        }

        return( sce )
    }
    else
        stop("Require package monocle to be installed to use this function.")
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

#' Merge SCESet objects
#'
#' Merge two SCESet objects that have the same features but contain different cells/samples.
#'
#' @param x an \code{\link{SCESet}} object
#' @param y an \code{\link{SCESet}} object
#' @param fdata_cols a character vector indicating which columns of featureData
#' of \code{x} and \code{y} should be retained. Alternatively, an integer or 
#' logical vector can be supplied to subset the column names of \code{fData(x)},
#' such that the subsetted character vector contains the columns to be retained.
#' Defaults to all shared columns between \code{fData(x)} and \code{fData(y)}.
#' @param pdata_cols a character vector indicating which columns of phenoData
#' of \code{x} and \code{y} should be retained. Alternatively, an integer or
#' logical vector to subset the column names of \code{pData(x)}. Defaults to
#' all shared columns between \code{pData(x)} and \code{pData(y)}.
#'
#' @details Existing cell-cell pairwise distances and feature-feature pairwise distances will not be valid for a merged \code{SCESet} object.
#' These entries are subsequently set to \code{NULL} in the returned object. 
#' Similarly, new \code{experimentData} will need to be added to the merged object.
#' 
#' If \code{fdata_cols} does not include the definition of feature controls, the control sets may not be defined in the output object.
#' In such cases, a warning is issued and the undefined control sets are removed from the \code{featureControlInfo} of the merged object.
#'
#' It is also \emph{strongly} recommended to recompute all size factors using the merged object, and re-run \code{\link{normalize}} before using \code{exprs}.
#' For arbitrary \code{x} and \code{y}, there is no guarantee that the size factors (and thus \code{exprs}) are comparable across objects.
#'
#' @return a merged \code{SCESet} object combining data and metadata from \code{x} and \code{y}
#'
#' @export
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' mergeSCESet(example_sceset[, 1:20], example_sceset[, 21:40])
#'
#' ## with specification of columns of fData
#' example_sceset <- calculateQCMetrics(example_sceset)
#' mergeSCESet(example_sceset[, 1:20], example_sceset[, 21:40], fdata_cols = c(1, 7))
#'
#' ## with specification of columns of pData
#' mergeSCESet(example_sceset[, 1:20], example_sceset[, 21:40], pdata_cols = 1:6)
#' mergeSCESet(example_sceset[, 1:20], example_sceset[, 40], pdata_cols = 3)
#'
mergeSCESet <- function(x, y, fdata_cols = NULL, pdata_cols = NULL) {
    if (!is(x,'SCESet')) stop('x must be of type SCESet')
    if (!is(y,'SCESet')) stop('y must be of type SCESet')
    if (!identical(featureNames(x), featureNames(y))) stop("feature names of x and y must be identical")

    for (sl in c("lowerDetectionLimit", "logExprsOffset", "featureControlInfo")) { 
        if (!identical(slot(x, sl), slot(y, sl))) 
            stop(sprintf("x and y do not have the same %s", sl))
    }

    ## check consistent fData
    if (is.null(fdata_cols)) { 
        fdata_cols <- intersect(colnames(fData(x)), colnames(fData(y)))
    } else if (!is.character(fdata_cols)) {
        fdata_cols <- colnames(fData(x))[fdata_cols]
    }
    fdata1 <- fData(x)[, fdata_cols, drop = FALSE]
    fdata2 <- fData(y)[, fdata_cols, drop = FALSE]
    if (!identical(fdata1, fdata2))
        stop("specified featureData columns are not identical for x and y")
    new_fdata <- as(fdata1, "AnnotatedDataFrame")

    ## combine pData
    if (is.null(pdata_cols)) { 
        pdata_cols <- intersect(colnames(pData(x)), colnames(pData(y)))
    } else if (!is.character(pdata_cols)) {
        pdata_cols <- colnames(pData(x))[pdata_cols]
    }
    pdata_x <- pData(x)[, pdata_cols, drop = FALSE]
    pdata_y <- pData(y)[, pdata_cols, drop = FALSE]
    if (!identical(colnames(pdata_x), colnames(pdata_y))) 
        stop("phenoData column names are not identical for x and y")
    if (ncol(pdata_x)) { 
        new_pdata <- rbind(pdata_x, pdata_y)
    } else {
        new_pdata <- data.frame(row.names = c(rownames(pdata_x), 
                                              rownames(pdata_y)))
    }
    new_pdata <- as(new_pdata, "AnnotatedDataFrame")

    ## combine exprsData
    new_exprs <- Biobase::combine(exprs(x), exprs(y))

    ## new SCESet
    merged_sceset <- newSCESet(exprsData = new_exprs, featureData = new_fdata,
                               phenoData = new_pdata,
                               lowerDetectionLimit = x@lowerDetectionLimit,
                               logExprsOffset = x@logExprsOffset)

    ## checking that the controls actually exist in the merged object
    all.fnames <- .fcontrol_names(x)
    discard <- logical(length(all.fnames))
    for (f in seq_along(all.fnames)) {
        fc <- all.fnames[f]
        which.current <- fData(merged_sceset)[[paste0("is_feature_control_", fc)]]
        if (is.null(which.current)) {
            warning(sprintf("removing undefined feature control set '%s'", fc))
            discard[f] <- TRUE
        }
    }
    featureControlInfo(merged_sceset) <- featureControlInfo(x)[!discard,]

    ## add remaining assayData to merged SCESet
    assay_names <- intersect(names(Biobase::assayData(x)),
                             names(Biobase::assayData(y)))
    for (assaydat in assay_names) {
        new_dat <- Biobase::combine(get_exprs(x, assaydat), get_exprs(y, assaydat))
        set_exprs(merged_sceset, assaydat) <- new_dat
    }
    merged_sceset
}



################################################################################
## writeSCESet

#' Write an SCESet object to an HDF5 file
#'
#' @param object \code{\link{SCESet}} object to be writted to file
#' @param file_path path to written file containing data from SCESet object
#' @param type character string indicating type of output file. Default is "HDF5".
#' @param overwrite_existing logical, if a file of the same name already exists
#' should it be overwritten? Default is \code{FALSE}.
#'
#' @details Currently writing to HDF5 files is supported. The \pkg{\link{rhdf5}}
#' package is used to write data to file and can be used to read data from HDF5
#' files into R. For further details about the HDF5 data format see
#' \url{https://support.hdfgroup.org/HDF5/}.
#'
#' @return Return is \code{NULL}, having written the \code{SCESet} object to file.
#' @export
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#'
#' \dontrun{
#' writeSCESet(example_sceset, "test.h5")
#' file.remove("test.h5")
#' }
#'
writeSCESet <- function(object, file_path, type = "HDF5", overwrite_existing = FALSE) {
    if (!is(object,'SCESet')) stop('object must be of type SCESet')
    if (file.exists(file_path) && !overwrite_existing)
        stop("To overwrite an existing file use argument overwrite_existing=TRUE")
    if (file.exists(file_path))
        file.remove(file_path)
    if (type == "HDF5") {
        rhdf5::H5close()
        rhdf5::h5createFile(file_path)
        tryCatch({
            rhdf5::h5write(featureNames(object), file = file_path,
                           name = "featureNames")
            rhdf5::h5write(cellNames(object), file = file_path, name = "cellNames")
            rhdf5::h5write(object@logExprsOffset, file = file_path,
                           name = "logExprsOffset")
            rhdf5::h5write(object@lowerDetectionLimit, file = file_path,
                           name = "lowerDetectionLimit")
            if (ncol(pData(object)) > 0)
                rhdf5::h5write(pData(object), file = file_path, name = "phenoData",
                               write.attributes = FALSE)
            if (ncol(fData(object)) > 0)
                rhdf5::h5write(fData(object), file = file_path, name = "featureData",
                               write.attributes = FALSE)
            rhdf5::h5createGroup(file_path, "assayData")
            for (assay in names(Biobase::assayData(object))) {
                group_set <- paste0("assayData/", assay)
                rhdf5::h5write(get_exprs(object, assay), file = file_path,
                               name = group_set,
                               write.attributes = FALSE)
            }
            rhdf5::H5close()
        }, finally = rhdf5::H5close())
    } else
        stop("HDF5 is the only format currently supported. See also saveRDS() to write to an object readable with R.")
}
