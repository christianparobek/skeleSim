setClassUnion("logOrNULL", c("logical", "NULL"))
setClassUnion("funcOrNULL", c("function", "NULL"))
setClassUnion("matrOrNULL", c("matrix", "NULL"))
setClassUnion("dfOrNULL", c("data.frame", "NULL"))

#' @name fastsimcoal.classes
#' @title fastsimcoal Parameters Class
#' @description An S4 class storing parameters specific to fastsimcoal
#'
#' @slot fastsimcoal.exec character string for the fastsimcoal command line
#'   executable.
#' @slot sample.times a vector giving the number of generations in the past
#'   at which samples are taken.
#' @slot growth.rate a vector giving the growth rate of each population.
#' @slot hist.ev a matrix describing historical events.
#' @slot num.chrom a value giving the number of chromosomes that the
#'   \code{locus.params} marker specifications should be copied for. If
#'   \code{NULL}, then chromosome assignment is taken from the
#'   \code{chromosome} column in \code{locus.params}. Any non-\code{NULL}
#'   integer will cause the \code{chromosome} column to be ignored.
#' @slot locus.params a data.frame giving the parameters for each locus.
#' @slot inf.site.model a logical indicating whether to use the infinite
#'   site model in fastsimcoal.
#'
setClass(
  Class = "fastsimcoal.params",
  slots = c(
    fastsimcoal.exec = "character", sample.times = "intOrNum",
    growth.rate = "intOrNum", hist.ev = "matrOrNULL", num.chrom = "intOrNum",
    locus.params = "dfOrNULL", inf.site.model = "logOrNULL"
  ),
  prototype = c(
    fastsimcoal.exec = "fsc252", sample.times = NULL, growth.rate = NULL,
    hist.ev = NULL, num.chrom = NULL, locus.params = NULL, inf.site.model = NULL
  ),
  validity = function(object) {
    return(TRUE)
  }
)

setMethod("initialize", "fastsimcoal.params",
  function(.Object, fastsimcoal.exec = "fsc252", sample.times = NULL,
           growth.rate = NULL, hist.ev = NULL,
           locus.params = NULL, num.chrom = NULL, inf.site.model = NULL
  ) {
    .Object@fastsimcoal.exec <- fastsimcoal.exec
    .Object@sample.times <- sample.times
    .Object@growth.rate <- growth.rate
    .Object@hist.ev <- hist.ev
    .Object@locus.params <- locus.params
    .Object@num.chrom <- num.chrom
    .Object@inf.site.model <- inf.site.model
    .Object
  }
)

#' @rdname fastsimcoal.classes
#' @export
#'
setGeneric("loadLocusParams<-", function(x, value) standardGeneric("loadLocusParams<-"))

#' @rdname fastsimcoal.classes
#' @importFrom methods validObject
#' @export
#'
setMethod("loadLocusParams<-", "fastsimcoal.params", function(x, value) {
  if(!inherit(value, "fsc.locusParams")) {
    stop("'value' must be a one row data.frame from either 'fsc.locus.dna', 'fsc.locus.snp', or 'fsc.locus.msat'")
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