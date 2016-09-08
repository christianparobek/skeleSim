#' @title Check parameters for rmetasim
#' @description Check parameters for rmetasim
#'
#' @param params a \linkS4class{skeleSim.params} object.
#'
#' @export
#'
rms.scenarioCheck <- function(params) {
  # check that sample times and growth rates are of length number of populations, and
  #   that historical events matrix converges
  results <- sapply(params@scenarios, function(sc) {

    c(
      dem.matr.same.dims = all(dim(sc@simulator.params@surv.matr)==
                                 dim(sc@simulator.params@repr.matr)),
      freqs.leng.num.alleles = length(sc@simulator.params@allele.freqs)==
                                  length(sc@simulator.params@num.alleles),
      nall.leng.num.loci = length(sc@simulator.params@num.alleles)==sc@num.loci,
      afrqs.leng.num.loci = length(sc@simulator.params@allele.freqs)==sc@num.loci,
      mut.leng.num.loci = length((sc@mut.rate))==sc@num.loci
    )
  })
  return(results)
}
