#' @name fsc.run
#' @title Run fastsimcoal
#' @description Run fastsimcoal
#'
#' @param params a \linkS4class{skeleSim.params} object.
#'
#' @return a modified \linkS4class{skeleSim.params} object with the results of
#'   a fastsimcoal run.
#'
#' @export
#'
fsc.run <- function(params) {
  label <- currentLabel(params)
  sc <- currentScenario(params)

  # Check that folder is empty
  if(file.exists(params@wd)) for(f in dir(label, full.names = T)) file.remove(f)

  params@rep.sample <- fastsimcoal(
    pop.info = sc@simulator.params@pop.info,
    locus.params = sc@simulator.params@locus.params,
    mig.rates = sc@migration,
    hist.ev = sc@simulator.params@hist.ev,
    label = label,
    quiet = params@quiet,
    exec = sc@simulator.params@fastsimcoal.exec,
    num.cores = 1
  )

  return(params)
}