# A .Call function with inbuilt checking for a returned error message.

.checkedCall <- function(cxxfun, ...) {
    out <- .Call(cxxfun, ...)
    if (is.character(out)) { stop(out) }
    return(out)
}
