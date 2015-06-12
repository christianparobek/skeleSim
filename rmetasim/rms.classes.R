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
#' @slot num.stg a single integer for number of stages
#' @slot self.rate a single numeric for selfing rate
#' @slot surv.matr a matrix of numerics for survival transitions between stages
#' @slot repr.matr a matrix of integers for reproductive (females)
#' @slot male.matr a matrix of integers for male contribution (males)
#' @slot carrying a single integer for carrying capacity
#' @slot init.pop.sizes a vector of integers for initial starting population census sizes (length number of pop'ns)
#' @slot num.alleles a vector of integers for number of alleles per locus (length number of loci)
#' @slot allele.freq a vector of numerics for initial allele frequencies ???
#' @slot num.gen a single integer of number of generations to run
#'
setClass(
  Class = "rmetasim.params",
  slots = c(num.stg = "intOrNum", self.rate = "intOrNum",
            surv.matr = "intOrNum", repr.matr = "intOrNum",
            male.matr = "intOrNum", carrying = "intOrNum",
            init.pop.sizes = "intOrNum", num.alleles = "intOrNum",
            
            num.gen = "intOrNum"
  ),
  prototype = c(num.stg = NULL, self.rate = NULL,
            surv.matr = NULL, repr.matr = NULL,
            male.matr = NULL, carrying = NULL,
            init.pop.sizes = NULL, num.alleles = NULL,
            
            num.gen = NULL
  )
)
