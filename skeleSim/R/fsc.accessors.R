#' @name fsc.accessors
#' @title fastsimcoal accessor functions
#' @description fastsimcoal accessor functions
#'
NULL

#' @rdname fsc.locusParams
#' @param x a \linkS4class{fastsimcoal.params} object.
#' @export
#'
setGeneric("loadLocusParams<-", function(x, value) standardGeneric("loadLocusParams<-"))


#' @rdname fsc.locusParams
#' @param value locus parameter data.frame to load.
#' @importFrom methods validObject
#' @export
#'
setMethod("loadLocusParams<-", "fastsimcoal.params", function(x, value) {
  if(!inherits(value, "fsc.locusParams")) {
    stop("'value' must be a data.frame from either 'fsc.locus.dna', 'fsc.locus.snp', or 'fsc.locus.msat'")
  }

  if(is.null(x@locus.params)) {
    x@locus.params <- value
  } else {
    x@locus.params <- rbind(x@locus.params, value)
  }
  x@locus.params <- x@locus.params[order(x@locus.params$chromosome), ]
  validObject(x)
  x
})