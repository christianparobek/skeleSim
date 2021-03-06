#' @title Plot all results
#' @description Plot all results
#'
#' @param params a \linkS4class{skeleSim.params} object.
#'
#' @importFrom reshape2 melt
#' @import ggplot2
#'
plot.all.stats <- function(params) {
  results.df <- as.data.frame(params@analysis.results)
  num.sc <- length(params@scenarios)
  names.stats <- colnames(params@analysis.results)
  colnames(results.df) <- c(names.stats)
  num.stats <- length(names.stats) - 1

  #plotting option 1
  results.df <- melt(
    results.df, id = "scenario", measure.vars = c(names.stats[-1])
  )
  ggplot(results.df, aes_string("value")) +
    geom_density() +
    facet_wrap(scenario ~ variable, ncol = num.stats, scales = "free")

  #histogram plot
  ggplot(results.df, aes_string("value")) +
    geom_histogram() +
    facet_grid(variable ~ scenario, scales = "free")
 }
