#' @title Current Scenario
#' @description Return the parameters for the current scenario.
#'
#' @param param a \linkS4class{skeleSim.params} object.

currentScenario <- function(params) {
  params@scenarios[[params@current.scenario]]
}
