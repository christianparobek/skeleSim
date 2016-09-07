#' @title Run rmetasim
#' @description Run rmetasim
#'
#' @param params a \linkS4class{skeleSim.params} object.
#'
#' @return a modified \linkS4class{skeleSim.params} object with the results of
#'   a rmetasim run.
#'
#' @importFrom rmetasim is.landscape landscape.simulate landscape.sample landscape.make.genind
#' @export
#'
rms.run <- function(params) {
  sc <- currentScenario(params)
  #create rmetasim landscape skeletonland
  skeleland<-rms.init.landscape(
    num.pops = sc@num.pops,
    carrying = sc@pop.size,
    mig.rates = sc@migration[[1]],
    num.loc = sc@num.loci,
    loc.type = sc@locus.type,
    mut.rate = sc@mut.rate,
    seq.length = sc@sequence.length,
    num.stgs = sc@simulator.params@num.stgs,
    selfing = sc@simulator.params@selfing,
    surv.matr = sc@simulator.params@surv.matr,
    repr.matr = sc@simulator.params@repr.matr,
    male.matr = sc@simulator.params@male.matr,
    init.pop.sizes = round(
      sum(sc@pop.size) * rep(1/(sc@num.pops*sc@simulator.params@num.stgs),
                             sc@num.pops*sc@simulator.params@num.stgs)
    ),
    num.alleles = sc@simulator.params@num.alleles,
    allele.freqs = sc@simulator.params@allele.freqs)

  #check is landscape ok
  if (!is.landscape(skeleland)) stop("landscape not cool")

  print(table(skeleland$individuals[,1]))
  print(paste("simulate for ",sc@simulator.params@num.gen,"years"))

  #run a number of generations
  skeleland<-landscape.simulate(skeleland, sc@simulator.params@num.gen)

  print("was able to simulate a rmetasim landscape")

  #take samples
  skeleland_samp<-landscape.sample(skeleland, ns=24)  ###need to improve
  #print(skeleland_samp)

  #print("about to check sampling")
  #check is landscape ok still
  if (!is.landscape(skeleland_samp)) stop("landscape not cool after sampling")

  #print(str(skeleland_samp))

  #  print("was able to sample")

  #  print("bout to convert")
  #now store the results
  save(file="tmpskeleland.rda",skeleland_samp)
  params@rep.sample<- rms.convert(skeleland_samp, sc@locus.type)
  print("returned from convert; returning params obj")
  return(params)
}


#' @title Write metasim file
#' @description Write metasim landscape script to disk (one per scenario)
#'
#' @param params a \linkS4class{skeleSim.params} object.
#'
#' @return Nothing
#'
#' @importFrom rmetasim is.landscape landscape.simulate landscape.sample landscape.make.genind
#' @export
#'
rms.write <- function(params) {
  numsc <- length(params@scenarios)
  for (i in 1:numsc) {
    fn <- paste0(gsub(" ","",params@title),"-",i,"-landscape-fun.R")
    params@current.scenario <- i
    sc <- currentScenario(params)
    outfile <- rms.init.landscape.func(
      num.pops = sc@num.pops,
      carrying = sc@pop.size,
      mig.rates = sc@migration[[1]],
      num.loc = sc@num.loci,
      loc.type = sc@locus.type,
      mut.rate = sc@mut.rate,
      seq.length = sc@sequence.length,
      num.stgs = sc@simulator.params@num.stgs,
      selfing = sc@simulator.params@selfing,
      surv.matr = sc@simulator.params@surv.matr,
      repr.matr = sc@simulator.params@repr.matr,
      male.matr = sc@simulator.params@male.matr,
      init.pop.sizes = round(
        sum(sc@pop.size) * rep(1/(sc@num.pops*sc@simulator.params@num.stgs),
                               sc@num.pops*sc@simulator.params@num.stgs)
      ),
      num.alleles = sc@simulator.params@num.alleles,
      allele.freqs = sc@simulator.params@allele.freqs
    )
    cat(file=fn,outfile)
  }
}
