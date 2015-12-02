#' @title Check general simulation parameters
#' @description Check general simulation parameters
#'
#' @param params a \linkS4class{skeleSim.params} object
#'
#' @export
#'
gen.scenario.check <- function(params) {
  #check that number of populations is same as length of pop sizes
  #check that number of populations is same as length of sample sizes
  #check that migration matrix is square with sises equal to number pops
  #check that num.pops has to be 1 or greater
  #check that num.loci has to be 1 or greater
  #check that mut.rate between 0 (inclusive) and 1 (exclusive)
  #check that at least one sample sizes has to be greater than 0
  #check that migration matrix all between 0 and 1
  #check that all migration matrix diagonals are 0

  results.check <- sapply(params@scenarios, function(sc) {
    #This will make a vector of TRUE/ FALSE
    c(nsizes.eq.npops = length(sc@pop.size) == sc@num.pops,
      nsamps.eq.npops = length(sc@sample.size) == sc@num.pops,
      is.mig.square = sapply(sc@migration, function(mig) {
        nrow(mig) == ncol(mig) & nrow(mig) == sc@num.pops
        }),
      at.lst.1.pop = sc@num.pops >= 1,
      at.lst.1.loc = sc@num.loci >= 1,
      mut.rate.ok = all((sc@mut.rate>=0)&(sc@mut.rate<1)),
      at.lst.1.samp = min(sc@sample.size)>0,
      mig.bet.0.1 = sapply(sc@migration, function(mig) {
        all((mig>=0)&(mig<=1))
      }),
      mig.diag.eq.0 = sapply(sc@migration, function(mig) {
        all(diag(mig)==0)
      })
    )
  })
  return(results.check)
}
