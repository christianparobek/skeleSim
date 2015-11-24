#' @title Current scenario
#'
#' @param param a \linkS4class{skeleSim.params} object.
#'
#' @return the parameters for the current scenario.
#'
currentScenario <- function(params) {
  params@scenarios[[params@current.scenario]]
}
