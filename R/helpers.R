## A suite of useful, miscellaneous functions

#' #' Trim whitespace from start and end of a (vector of) string(s)
#'
#' @param x vector of strings (of length one or greater)
#' @return vector of strings with whitespace trimmed from start and end of each string
#' @export
#' @examples
#' trim("   This is a string   ")
trim <- function(x) gsub("^\\s+|\\s+$", "", x)


















