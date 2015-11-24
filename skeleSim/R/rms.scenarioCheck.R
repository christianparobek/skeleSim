rms.scenarioCheck <-
function(params) {
  # check that sample times and growth rates are of length number of populations, and
  #   that historical events matrix converges 
  results <- sapply(params@scenarios, function(sc) {
    
    c(
      dem.matr.same.dims = all(dim(sc@simulator.params@surv.matr)==
                                 dim(sc@simulator.params@repr.matr)),
      freqs.leng.num.alleles = length(sc@simulator.params@allele.freqs)==
                                  sc@simulator.params@num.alleles,
      mut.leng.num.loci = length((sc@mut.rate))==sc@num.loci
    )
  })
  return(results)
}
