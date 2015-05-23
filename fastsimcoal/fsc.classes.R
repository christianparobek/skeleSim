setClassUnion("logOrNULL", c("logical", "NULL"))
setClassUnion("listOrNULL", c("list","NULL"))
setClassUnion("charOrNULL", c("character", "NULL"))
setClassUnion("intOrNum", c("integer","numeric", "NULL"))
setClassUnion("funcOrNULL", c("function", "NULL"))
setClassUnion("posixOrNULL", c("POSIXct", "POSIXlt", "NULL"))

#' @title fastsimcoal Parameters Class
#' @description An S4 class storing parameters specific to fastsimcoal
#'
#' @rdname fastsimcoal.classes
#'
#' @slot sample.times a vector giving the number of generations in the past
#'   at which samples are taken
#' @slot growth.rate a vector giving the growth rate of each population
#' @slot hist.ev a matrix describing historical events
#' @slot locus.params a list giving the parameters for each locus.
#'
setClass(
  Class = "fastsimcoal.params",
  slots = c(sample.times = "intOrNum", growth.rate = "intOrNum",
            hist.ev = "intOrNum", locus.params = "intOrNum",
            inf.site.model = "logOrNULL"
  ),
  prototype = c(sample.times = NULL, growth.rate = NULL, hist.ev = NULL,
                locus.params = NULL, inf.site.model = NULL
  )
)