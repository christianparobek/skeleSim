#' @title Current label
#'
#' @param param a \linkS4class{skeleSim.params} object.
#'
#' @return a character label representing current scenario and replicate.
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
