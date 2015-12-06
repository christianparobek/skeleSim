#' @title GUI for skeleSim system
#' @description This function starts the shiny simulation control panel
#'
#'
#' @return NULL
#' @import shiny
#' @importFrom igraph graph.adjacency plot.igraph
#' @export
#' 
skeleSimGUI <- function() {
  appDir <- system.file("shiny", "skeleSimShinyGUI", package = "skeleSim")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `skeleSim`.", call. = FALSE)
  }
  runApp(appDir, display.mode = "normal")
}
