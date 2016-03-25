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

  locus.type <- c("dna","msat","snp")[which(c("DNA","MICROSAT")==sc@simulator.params@locus.params[1,1])]
  
  params@rep.sample <- fastsimcoal(
    pop.info = fscPopInfo(pop.size=sc@pop.size,
                          sample.size=sc@sample.size,
                          sample.times=sc@simulator.params@sample.times,
                          growth.rate=sc@simulator.params@growth.rate),
#    locus.params = sc@simulator.params@locus.params,
    locus.params = fscLocusParams(locus.type=locus.type,num.loci=1,mut.rate=sc@simulator.params@locus.params[,4],
                                  chromosome=1:dim(sc@simulator.params@locus.params)[1],
                                  sequence.length=ifelse(locus.type=="dna",sc@sequence.length,NULL)),
    mig.rates = sc@migration,
    hist.ev = sc@simulator.params@hist.ev,
    label = label,
    quiet = params@quiet,
    exec = sc@simulator.params@fastsimcoal.exec,
    num.cores = 1
  )

  return(params)
}
