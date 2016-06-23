# A .Call function with inbuilt checking for a returned error message.

#' @useDynLib scater, .registration=TRUE, .fixes="cxx_"
.checkedCall <- function(cxx_fun, ...) {
    out <- .Call(cxx_fun, ...)
    if (is.character(out)) { stop(out) }
    return(out)
}
