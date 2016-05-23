setClassUnion("matrOrNULL", c("matrix", "NULL"))

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
#'
#' @export
#'
fastsimcoal.params <- setClass(
  Class = "fastsimcoal.params",
  slots = c(
    fastsimcoal.exec = "character", sample.times = "intOrNum",
    growth.rate = "intOrNum", hist.ev = "matrOrNULL", num.chrom = "intOrNum",
    locus.params = "dfOrNULL"
  ),
  prototype = list(
    fastsimcoal.exec = "fsc252", sample.times = NULL, growth.rate = NULL,
    hist.ev = NULL, num.chrom = NULL, locus.params = NULL
  ),
  validity = function(object) {
    return(TRUE)
  }
)