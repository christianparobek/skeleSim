#' @title Current Label
#' @description Return a character label for the current scenario and replicate
#'
#' @param param a \linkS4class{skeleSim.params} object.

currentLabel <- function(params) {
  label <- paste(
    params@title,
    params@current.scenario,
    params@current.replicate,
    sep = "."
  )
  gsub("[[:punct:]]", ".", label)
}
