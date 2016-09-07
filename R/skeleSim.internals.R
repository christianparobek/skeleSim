#' @name skeleSim.internals
#' @title SkeleSim internal functions
#' @description SkeleSim internal functions
#'
#' @param params a \linkS4class{skeleSim.params} object.
#'
#' @return
#'   \tabular{ll}{
#'     \code{currentScenario} \tab the parameters for the current scenario.\cr
#'     \code{currentLabel} \tab a character label representing current scenario
#'       and replicate.\cr
#'   }
#'
#' @keywords internals
#'
currentScenario <- function(params) {
  params@scenarios[[params@current.scenario]]
}

#' @rdname skeleSim.internals
#' @keywords internals
#'
currentLabel <- function(params) {
  label <- paste(
    params@title,
    params@current.scenario,
    params@current.replicate,
    sep = "."
  )
  gsub("[[:punct:]]", ".", label)
}
