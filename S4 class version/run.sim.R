#' @title Run Simulation
#' @description Run simulator using specified parameters.
#'
#' @param params a complete \linkS4class{skeleSim.params} object.
#'
#' @note Function first tests validity of \code{params} and returns
#'   \code{NULL} if invalid.
#'
run.sim <- function(params) {
  if(!validObject(params, test = TRUE, complete = TRUE)) invisible(NULL)
  params <- sim.iterator(params)
  save(params, file = paste(params@title, ".last.rep.rdata", sep = ""))
  plot_all_stats(params)
  params <- summary_table_stats(params)
  cat("/n")
  print(params@summary.results)
  save(params, file = paste(params@title, ".results.rdata", sep = ""))
  invisible(params)
}