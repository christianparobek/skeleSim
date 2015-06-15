rms.run<-function(params){
  sc <- currentScenario(params)
  
  #create rmetasim landscape skeletonland
  skeleland<-rms.init.landscape(
      num.pops = sc@num.pops, carrying = sc@pop.size, sample.size = sc@sample.size, mig.rates = sc@migration, 
      num.loc = sc@num.loci, loc.type = sc@lucus.type, mut.rate = sc@mut.rate, seq.length = sc@sequence.length, 
      num.stgs = sc@simulator.params@num.stgs, selfing = sc@simulator.params@selfing, 
      surv.matr = sc@simulator.params@surv.matr, repr.matr = sc@simulator.params@repr.matr, 
      male.matr = sc@simulator.params@male.matr,
      init.pop.sizes = sc@simulator.params@init.pop.sizes, num.gen = sc@simulator.params@num.gen, 
      num.alleles = sc@simulator.params@num.alleles, allele.freqs = sc@simulator.params@allele.freqs)
  
  #run a number of generations
  skeleland<-landscape.simulate(skeleland, num.gen)

  #take samples
  skeleland_samp<-landscape.sample(skeleland, sample.size)
  
  #now store the results
  params@rep.sample<- rms.convert(skeleland_samp)

}
