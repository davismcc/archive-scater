rechunk <- function(incoming, redim, outfile=NULL, outname=NULL, byrow=TRUE) {
    # Setting names, if not present.
    if (is.null(outfile)) { 
        outfile <- HDF5Array::getHDF5DumpFile(for.use=TRUE)        
    } else {
        if (!h5createFile(outfile)) {
            stop("failed to create output HDF5 file")
        }
    }
    if (is.null(outname)) { 
        outname <- HDF5Array::getHDF5DumpName(for.use=TRUE)        
    }

    # Repacking the file.
    if (byrow) {
        chunk.dims <- c(1L, min(ncol(incoming), redim))
    } else {
        chunk.dims <- c(min(nrow(incoming), redim), 1L)
    }
    data.type <- type(incoming)
    .Call(cxx_rechunk_matrix, incoming@seed@file, incoming@seed@name, data.type,
          outfile, outname, redim, byrow, PACKAGE="scater")

    # Generating output. 
    HDF5Array::appendDatasetCreationToHDF5DumpLog(outfile, outname, dim(incoming), 
                                                  data.type, "?", chunk.dims)
    HDF5Array::HDF5Array(outfile, outname, data.type)
}

