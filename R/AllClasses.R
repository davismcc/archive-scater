### all classes defined for the scater package

################################################################################
### defining the SCESet class

#' The "Single Cell Expression Set" (SCESet)  class
#'
#' S4 class and the main class used by scater to hold single cell expression 
#' data. SCESet extends the basic Bioconductor ExpressionSet class.
#'
#' This class is initialized from a matrix of expression values.
#' 
#' Methods that operate on SCESet objects constitute the basic scater workflow.
#'
#' Thanks to the Monocle package for their CellDataSet class, which provided the 
#' inspiration and template for SCESet.
#'
#'@section Slots:
#'  \describe{
#'    \item{\code{logged}:}{Scalar of class \code{"logical"}, indicating whether 
#'    or not the expression data in the `exprs` slot have been log2-transformed
#'    or not.}
#'    \item{\code{lowerDetectionLimit}:}{Scalar of class \code{"numeric"}, 
#'    giving the lower limit for an expression value to be classified as 
#'    "expressed".}
#'}
#' @name SCESet
#' @rdname SCESet
#' @aliases SCESet-class
#' @exportClass SCESet
setClass( "SCESet",
          contains = "ExpressionSet",
          slots = c(logged = "logical",
                    lowerDetectionLimit = "numeric"),
          prototype = prototype(new("VersionedBiobase",
                                    versions = c(classVersion("ExpressionSet"),
                                                 SCESet = "0.1.1")))
)


