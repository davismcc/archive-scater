# ### defining all generic methods
 
#' Replacement for phenoData
#' 
#' @name pData<-
#' @param x object containing phenoData to replace
#' @param value AnnotatedDataFrame to replace existing phenoData
#' @export
#' @docType methods
#' @rdname pData-replace
#' @importFrom Biobase pData<-
setGeneric("pData<-", signature = signature("x", "value"),
           function(x, value) {standardGeneric("pData<-")})

#' Replacement for featureData
#' 
#' @name fData<-
#' @param x object containing featureData to replace
#' @param value AnnotatedDataFrame to replace existing featureData
#' @export
#' @docType methods
#' @rdname fData-replace
#' @importFrom Biobase fData<-
setGeneric("fData<-", signature = signature("x", "value"),
           function(x, value) {standardGeneric("fData<-")})


#' @name isExpr
#' @export
#' @docType methods
#' @rdname isExpr 
setGeneric("isExpr", function(object) {standardGeneric("isExpr")})

#' @name isExpr<-
#' @export
#' @docType methods
#' @rdname isExpr
setGeneric("isExpr<-", function(object, value) {standardGeneric("isExpr<-")})


