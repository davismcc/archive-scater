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


#' @name is_exprs
#' @export
#' @docType methods
#' @rdname is_exprs
setGeneric("is_exprs", function(object) {standardGeneric("is_exprs")})

#' @name is_exprs<-
#' @export
#' @docType methods
#' @rdname is_exprs
setGeneric("is_exprs<-", function(object, value) {standardGeneric("is_exprs<-")})


#' @name norm_exprs
#' @export
#' @docType methods
#' @rdname norm_exprs
setGeneric("norm_exprs", function(object) {standardGeneric("norm_exprs")})

#' @name norm_exprs<-
#' @export
#' @docType methods
#' @rdname norm_exprs
setGeneric("norm_exprs<-", function(object, value) {standardGeneric("norm_exprs<-")})

#' @name stand_exprs
#' @export
#' @docType methods
#' @rdname stand_exprs
setGeneric("stand_exprs", function(object) {standardGeneric("stand_exprs")})

#' @name stand_exprs<-
#' @export
#' @docType methods
#' @rdname stand_exprs
setGeneric("stand_exprs<-", function(object, value) {standardGeneric("stand_exprs<-")})


#' @name tpm
#' @export
#' @docType methods
#' @rdname tpm
setGeneric("tpm", function(object) {standardGeneric("tpm")})

#' @name tpm<-
#' @export
#' @docType methods
#' @rdname tpm
setGeneric("tpm<-", function(object, value) {standardGeneric("tpm<-")})

#' @name cpm
#' @export
#' @docType methods
#' @rdname cpm
setGeneric("cpm", function(object) {standardGeneric("cpm")})

#' @name cpm<-
#' @export
#' @docType methods
#' @rdname cpm
setGeneric("cpm<-", function(object, value) {standardGeneric("cpm<-")})

#' @name fpkm
#' @export
#' @docType methods
#' @rdname fpkm
setGeneric("fpkm", function(object) {standardGeneric("fpkm")})

#' @name fpkm<-
#' @export
#' @docType methods
#' @rdname fpkm
setGeneric("fpkm<-", function(object, value) {standardGeneric("fpkm<-")})

#' @name bootstraps
#' @export
#' @docType methods
#' @rdname bootstraps
setGeneric("bootstraps", function(object) {standardGeneric("bootstraps")})

#' @name bootstraps<-
#' @export
#' @docType methods
#' @rdname bootstraps
setGeneric("bootstraps<-", function(object, value) {standardGeneric("bootstraps<-")})

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


#' @name reducedDimension
#' @export
#' @docType methods
#' @rdname reducedDimension
setGeneric("reducedDimension", function(object) {
    standardGeneric("reducedDimension")
})

#' @name reducedDimension<-
#' @export
#' @docType methods
#' @rdname reducedDimension
setGeneric("reducedDimension<-", function(object, value) {
    standardGeneric("reducedDimension<-")
})


#' @name redDim
#' @export
#' @docType methods
#' @rdname reducedDimension
setGeneric("redDim", function(object) {
    standardGeneric("redDim")
})

#' @name redDim<-
#' @export
#' @docType methods
#' @rdname reducedDimension
setGeneric("redDim<-", function(object, value) {
    standardGeneric("redDim<-")
})

#' @name plotReducedDim
#' @export
#' @docType methods
#' @rdname plotReducedDim
setGeneric("plotReducedDim", function(object, ...) {
    standardGeneric("plotReducedDim")
})

#' @name plotExpression
#' @export
#' @docType methods
#' @rdname plotExpression
setGeneric("plotExpression", function(object, ...) {
    standardGeneric("plotExpression")
})


