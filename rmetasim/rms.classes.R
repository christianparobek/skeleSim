setClassUnion("logOrNULL", c("logical", "NULL"))
setClassUnion("listOrNULL", c("list","NULL"))
setClassUnion("charOrNULL", c("character", "NULL"))
setClassUnion("intOrNum", c("integer","numeric", "NULL"))
setClassUnion("funcOrNULL", c("function", "NULL"))
setClassUnion("posixOrNULL", c("POSIXct", "POSIXlt", "NULL"))

#' @title rmetasim Parameters Class
#' @description An S4 class storing parameters specific to rmetasim
#'
#' @rdname rmetasim.classes
#'
#' @slot num.stg number of stages
#' @slot self.rate selfing rate
#' @slot surv.matr survival matrix
#' @slot repr.matr reproductive matrix (females)
#' @slot male.matr male contribution matrix (males)
#' @slot carrying carrying capacity
#' @slot init.pop.sizes initial starting population census
#' @slot num.alleles number of alleles per locus
#' @slot allele.freq initial allele frequencies
#' @slot num.gen number of generations to run
#'
setClass(
  Class = "fastsimcoal.params",
  slots = c(num.stg = "intOrNum", self.rate = "intOrNum",
            surv.matr
  
  ),
  prototype = c(num.stg = NULL, self.rate = NULL,
  
  
  )
)
