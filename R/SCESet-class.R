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
#'    \item{\code{counts}:}{Matrix of class \code{"numeric"}, containing 
#'    the raw feature-level count data.}
#'    \item{\code{expr_vals}:}{Matrix of class \code{"numeric"}, containing 
#'    the feature-level expression values to be used for analyses. Typically 
#'    these will be count-per-million values derived from the count data, but 
#'    could be a different measure of expression. This is a default slot in the 
#'    Bioconductor ExpressionSet class and can be accessed with the \code{exprs}
#'    function.}
#'    \item{\code{cell_info}:}{Dataframe of containing cell metadata information
#'    used for QC and modelling.}  
#'  }
#'
#' @name SCESet
#' @rdname SCESet
#' @aliases SCESet-class
#' @exportClass SCESet
setClass( "SCESet",
          contains = "ExpressionSet",
          slots = c(counts = "matrix",
                    expr_vals = "matrix",
                    cell_info = "data.frame"),
          prototype = prototype(new("VersionedBiobase",
                                    versions = c(classVersion("ExpressionSet"),
                                                 SCESet = "0.1.0")))
)



