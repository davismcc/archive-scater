#' @name norm_exprs
#' @export
#' @docType methods
#' @return a matrix of normalised expression data
#' @rdname norm_exprs
setGeneric("norm_exprs", function(object) standardGeneric("norm_exprs"))

#' @name norm_exprs<-
#' @export
#' @docType methods
#' @rdname norm_exprs
setGeneric("norm_exprs<-", function(object, value) standardGeneric("norm_exprs<-"))

#' @name stand_exprs
#' @export
#' @docType methods
#' @return a matrix of standardised expressiond data
#' @rdname stand_exprs
setGeneric("stand_exprs", function(object) standardGeneric("stand_exprs"))

#' @name stand_exprs<-
#' @export
#' @docType methods
#' @rdname stand_exprs
setGeneric("stand_exprs<-", function(object, value) standardGeneric("stand_exprs<-"))

#' @name bootstraps
#' @export
#' @docType methods
#' @rdname bootstraps
setGeneric("bootstraps", function(object) standardGeneric("bootstraps"))

#' @name bootstraps<-
#' @export
#' @docType methods
#' @rdname bootstraps
setGeneric("bootstraps<-", function(object, value) standardGeneric("bootstraps<-"))

#' @name tpm
#' @export
#' @docType methods
#' @return a matrix of transcripts-per-million data
#' @rdname tpm
setGeneric("tpm", function(object) standardGeneric("tpm"))

#' @name tpm<-
#' @export
#' @docType methods
#' @rdname tpm
setGeneric("tpm<-", function(object, value) standardGeneric("tpm<-"))

#' @name cpm
#' @export
#' @docType methods
#' @return a matrix of counts-per-million values
#' @rdname cpm
setGeneric("cpm", function(object) standardGeneric("cpm"))

#' @name cpm<-
#' @export
#' @docType methods
#' @rdname cpm
setGeneric("cpm<-", function(object, value) standardGeneric("cpm<-"))

#' @name fpkm
#' @export
#' @docType methods
#' @return a matrix of FPKM values
#' @rdname fpkm
setGeneric("fpkm", function(object) standardGeneric("fpkm"))

#' @name fpkm<-
#' @export
#' @docType methods
#' @rdname fpkm
setGeneric("fpkm<-", function(object, value) standardGeneric("fpkm<-"))


