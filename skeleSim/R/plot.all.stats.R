#' @title Plot all results
#' @description Plot all results
#'
#' @param params a \linkS4class{skeleSim.params} object.
#'
#' @importFrom reshape2 melt
#' @import ggplot

plot.all.stats <- function(params) {
  results.datafr<-as.data.frame(params@analysis.results)
  num.sc <- length(params@scenarios)
  names.stats<-colnames(params@analysis.results)
  colnames(results.datafr)<-c(names.stats)
  num.stats<-length(names.stats)-1

  #plotting option 1
  results.melted<-melt(results.datafr,id="scenario",measure.vars=c(names.stats[-1]))
  ggplot(results.melted, aes(value)) +
      geom_density() + facet_wrap(scenario~variable,
      ncol = num.stats, scales = "free")

  #histogram plot
  ggplot(results.melted, aes(value)) + geom_histogram() +
  facet_grid(variable~scenario, scales="free")

 }
