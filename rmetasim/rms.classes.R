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
#' @slot num.gen number of generations to run
#'
setClass(
  Class = "fastsimcoal.params",
  slots = c(num.stg = "intOrNum", self.rate = "intOrNum",
  
  
  ),
  prototype = c(num.stg = NULL, self.rate = NULL,
  
  
  )
)
