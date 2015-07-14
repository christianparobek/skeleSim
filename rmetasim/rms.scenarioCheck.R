rms.scenarioCheck <- function(params) {
  # check that sample times and growth rates are of length number of populations, and
  #   that historical events matrix converges 
  results <- sapply(params@scenarios, function(sc) {
    
    c(
      dem.matr.same.dims <- dim(sc@simulator.params@surv.matr)==dim(sc@simulator.params@repr.matr)
    )
  })
  return(results)
}
