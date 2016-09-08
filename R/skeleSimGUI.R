#' @title GUI for skeleSim system
#' @description This function starts the shiny simulation control panel
#'
#' @param launch.browser If true, the system's default web browser will be
#'   launched automatically after the app is started.
#'
#' @return NULL
#'
#' # markdown and shinyFiles import added for shiny app
#' @import markdown shiny shinyFiles
#' @importFrom igraph graph.adjacency plot.igraph
#'
#' @export
#'
skeleSimGUI <- function(launch.browser = TRUE) {
  appDir <- system.file("shiny", "skeleSimShinyGUI", package = "skeleSim")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `skeleSim`.", call. = FALSE)
  }
  runApp(appDir, display.mode = "normal", launch.browser = launch.browser)
}