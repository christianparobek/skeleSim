fsc.run <- function(params) {
  label <- currentLabel(params)
  sc <- currentScenario(params)

  # Check that folder is empty
  if(file.exists(params@wd)) for(f in dir(label, full.names = T)) file.remove(f)

  # modify Ne and sample size if not sequence data to account for
  #   haploid -> diploid conversion of output
  size.mult <- if(sc@locus.type == "sequence") 1 else 2

  # Write fastsimcoal input file
  file <- fsc.write(
    num.pops = sc@num.pops,
    Ne = sc@pop.size * size.mult,
    sample.size = sc@sample.size * size.mult,
    sample.times = sc@simulator.params@sample.times,
    growth.rate = sc@simulator.params@growth.rate,
    mig.rates = sc@migration,
    num.chrom = sc@num.loci,
    hist.ev = sc@simulator.params@hist.ev,
    locus.params = sc@simulator.params@locus.params,
    label = label
  )

  # Run fastsimcoal
  err <- system(paste(
    sc@simulator.params@fastsimcoal.exec, "-i", file, "-n 1",
    ifelse(sc@simulator.params@inf.site.model, "-I", ""),
    ifelse(params@quiet, "-q", "")
  ), intern = F)

  if(err == 0) {
    if(!params@quiet) cat("fastsimcoal exited normally\n")
  } else {
    stop("fastsimcoal exited with error ", err, "\n")
  }

  arp.file <- file.path(label, paste(label, "_1_1.arp", sep = ""))
  params@rep.sample <- fsc.read(arp.file, params)
  params
}

fsc.scenarioCheck <- function(params) {
  # check that sample times and growth rates are of length number of populations, and
  #   that historical events matrix converges 
  results <- sapply(params@scenarios, function(sc) {
    
    c(
      samptime.eq.npops = length(sc@simulator.params@sample.times) == sc@num.pops,
      growth.eq.npops = length(sc@simulator.params@growth.rate) == sc@num.pops,
      hist.ev.good = fsc.histEvConverges(
                            sc@simulator.params@hist.ev,
                            sc@pop.size,
                            sc@simulator.params@growth.rate,
                            length(sc@migration))
    )
  })
  return(results)
}
