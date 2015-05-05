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


#' @name isExprs
#' @export
#' @docType methods
#' @rdname isExprs
setGeneric("isExprs", function(object) {standardGeneric("isExprs")})

#' @name isExprs<-
#' @export
#' @docType methods
#' @rdname isExprs
setGeneric("isExprs<-", function(object, value) {standardGeneric("isExprs<-")})


#' @name cellPairwiseDistances
#' @export
#' @docType methods
#' @rdname cellPairwiseDistances
setGeneric("cellPairwiseDistances", function(object) {
    standardGeneric("cellPairwiseDistances")
    })

#' @name cellPairwiseDistances<-
#' @export
#' @docType methods
#' @rdname cellPairwiseDistances
setGeneric("cellPairwiseDistances<-", function(object, value) {
    standardGeneric("cellPairwiseDistances<-")
    })


#' @name cellDist
#' @export
#' @docType methods
#' @rdname cellPairwiseDistances
setGeneric("cellDist", function(object) {
    standardGeneric("cellDist")
})

#' @name cellDist<-
#' @export
#' @docType methods
#' @rdname cellPairwiseDistances
setGeneric("cellDist<-", function(object, value) {
    standardGeneric("cellDist<-")
})


#' @name genePairwiseDistances
#' @export
#' @docType methods
#' @rdname genePairwiseDistances
setGeneric("genePairwiseDistances", function(object) {
    standardGeneric("genePairwiseDistances")
})

#' @name genePairwiseDistances<-
#' @export
#' @docType methods
#' @rdname genePairwiseDistances
setGeneric("genePairwiseDistances<-", function(object, value) {
    standardGeneric("genePairwiseDistances<-")
})


#' @name geneDist
#' @export
#' @docType methods
#' @rdname genePairwiseDistances
setGeneric("geneDist", function(object) {
    standardGeneric("geneDist")
})

#' @name geneDist<-
#' @export
#' @docType methods
#' @rdname genePairwiseDistances
setGeneric("geneDist<-", function(object, value) {
    standardGeneric("geneDist<-")
})



