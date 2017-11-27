#' @title Initialize a landscape object
#' @description Initialize a landscape object
#'
#' @param num.pops number of populations to simulate
#' @param carrying maximum population sizes for each population
#' @param sample.size size of sample to be pulled from each population
#' @param mig.rates a list of among-population migration matrices
#' @param num.loc number of independently segregating loci
#' @param loc.type sequence or microsatellite
#' @param mut.rate per gene mutation rate
#' @param seq.length if simulating a sequence, the length of the molecule
#' @param num.stgs number of demographic stages in a population
#' @param selfing selfing rate (must range from 0 [random mating] to 1 [complete selfing]
#' @param surv.matr within pop survival matrices
#' @param repr.matr within pop reproduction matrices
#' @param male.matr with pop male repro matrices
#' @param init.pop.sizes starting population sizes
#' @param num.gen number of generations to simulate
#' @param num.alleles vector of the number of alleles per locus
#' @param allele.freqs list of allele freqs for each locus (range 0-1)
#'
#' @importFrom rmetasim landscape.new.empty landscape.new.floatparam
#'   landscape.new.intparam landscape.new.switchparam landscape.new.local.demo
#'   landscape.mig.matrix landscape.new.epoch landscape.new.locus landscape.new.individuals
#' @export
rms.init.landscape <- function(
  num.pops = NULL, carrying = NULL, sample.size = NULL, mig.rates = NULL,
  num.loc = NULL, loc.type = NULL, mut.rate = NULL, seq.length = NULL,
  num.stgs = NULL, selfing = NULL, surv.matr = NULL, repr.matr = NULL,
  male.matr = NULL, init.pop.sizes = NULL, num.gen = NULL, num.alleles = NULL,
  allele.freqs = NULL) {

  if(exists("skeletonland")) rm(skeletonland)
  skeletonland<-landscape.new.empty()
  #define selfing rate
  skeletonland<-landscape.new.floatparam(skeletonland,s=selfing)
  #Hard coded in current generation, current epoch, max number generations, max number individuals
  skeletonland<-landscape.new.intparam(skeletonland,h=num.pops,s=num.stgs,cg=0,ce=0,totgen=100000,maxland=100000)
  #Hard coded in multiple paternity (yes) and density dependence (no) parameters
  skeletonland<-landscape.new.switchparam(skeletonland,re=0,rd=0,mp=1,dd=0)

  #local matrices, will give same demography to each local population
  for (i in 1:num.pops) {
    skeletonland<-landscape.new.local.demo(skeletonland,surv.matr, repr.matr, male.matr)
  }

  #cross habitat matrices
  #epoch_s_matr<-matrix(0,nrow=4, ncol=4)
  #epoch_r_matr<-
  #epoch_m_matr<-matrix(0,nrow=4, ncol=4)
  epoch_s_matr <- matrix(0,nrow=(num.pops*num.stgs), ncol=(num.pops*num.stgs))
  epoch_r_matr <- landscape.mig.matrix(h=num.pops,s=num.stgs,mig.model="custom",R.custom=mig.rates)$R
  epoch_m_matr <- epoch_s_matr


  #no extinction allowed, hard coded
  skeletonland <- landscape.new.epoch(
    skeletonland,epochprob=1, epoch_s_matr, epoch_r_matr, epoch_m_matr,
    startgen=0,extinct=NULL,carry=carrying
  )

  #LOCI
  #Note that for SNPs, numalleles should be 2, allelesize only used for sequences
  #type = 0 IAM, type = 1 SMM type = 2 DNA sequence
  #assumes biparental transmission (transmission = 0)
  rms.locus.type = NULL

  #print(loc.type)

  if (loc.type == "SNP") {rms.locus.type = 2; seq.length = rep(1,num.loc)}
  if (loc.type %in% c("microsat","MICROSAT","microsatellite")) rms.locus.type = 1
  if (loc.type == "sequence") rms.locus.type = 2
  for (l in 1:num.loc) {
    if (rms.locus.type==2) {
      if (loc.type=="sequence") { #a sequence not snps
        if (l==1) {#only one sequence locus possible and it creates a maternally inherited haploid locus
          skeletonland<-landscape.new.locus(
            skeletonland, type=2, ploidy=1, mutationrate=mut.rate[l],
            transmission=1, numalleles=num.alleles[l],
            frequencies=allele.freqs[[l]]
          ) #temp remove allele size
        }
      } else { #we are dealing with snps (diploid, biparental short seuqnces)
        skeletonland<-landscape.new.locus(
          skeletonland, type=2, ploidy=2, mutationrate=mut.rate[l],
          transmission=0, numalleles=num.alleles[l],
          frequencies=allele.freqs[[l]],
          allelesize=1
        )
      }

    } else {
      skeletonland<-landscape.new.locus(
        skeletonland, type=1, ploidy=2, mutationrate=mut.rate[l],
        transmission=0, numalleles=num.alleles[l],
        frequencies=allele.freqs[[l]]
      )
    }
  }
  ###assumes population initial sizes all defined nicely by user
  landscape.new.individuals(skeletonland,init.pop.sizes)
}