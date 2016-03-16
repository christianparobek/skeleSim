setClassUnion("matrOrNULL", c("matrix", "NULL"))
setClassUnion("logOrNULL", c("logical", "NULL"))
setClassUnion("intOrNULL", c("integer", "NULL"))

#' @title fastsimcoal Parameters Class
#' @description An S4 class storing parameters specific to fastsimcoal
#'
#' @slot fastsimcoal.exec character string for the fastsimcoal command line
#'   executable.
#' @slot pop.info matrix of population sampling information created by the
#'   \code{\link{fscPopInfo}} function.
#' @slot locus.params data.frame specifying loci to simulate created by the
#'   \code{\link{fscLocusParams}} function.
#' @slot hist.ev matrix of historical events created by the
#'   \code{\link{fscHistEv}} function.
#'
#' @export
#'
fastsimcoal.params <- setClass(
  Class = "fastsimcoal.params",
  slots = c(
    fastsimcoal.exec = "character", pop.info = "matrOrNULL",
    hist.ev = "matrOrNULL", locus.params = "dfOrNULL",
    inf.site.model = "logOrNULL", growth.rate = "intOrNULL",
    sample.times = "intOrNULL"
  ),
  prototype = list(
    fastsimcoal.exec = "fsc252", pop.info = NULL, hist.ev = NULL,
    locus.params = NULL, inf.site.model=T,growth.rate=0,sample.times=NULL
  ),
  validity = function(object) {
    return(TRUE)
  }
)
