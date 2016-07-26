#' @title Check parameters for fastsimcoal
#' @description Check parameters for fastsimcoal
#'
#' @param params a \linkS4class{skeleSim.params} object.
#'
#' @export
#'
fsc.scenarioCheck <- function(params) {
  # check that sample times and growth rates are of length number of populations, and
  #   that historical events matrix converges
  results <- sapply(params@scenarios, function(sc) {
    c(
      samptime.eq.npops = length(sc@simulator.params@sample.times) == sc@num.pops,
      growth.eq.npops = length(sc@simulator.params@growth.rate) == sc@num.pops,
      hist.ev.good = fsc.histEvCheck(
        hist.ev = sc@simulator.params@hist.ev,
        pop.size = sc@pop.size,
        growth.rate = sc@simulator.params@growth.rate,
        num.mig.mats = length(sc@migration)
      )
    )
  })
  return(results)
}
