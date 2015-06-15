# initiate an rmetasim landscape with parameters

rms.init.landscape <- function(num.pops = NULL, carrying = NULL, sample.size = NULL, mig.rates = NULL, 
                            num.loc = NULL, loc.type = NULL, mut.rate = NULL, seq.length = NULL, 
                            num.stgs = NULL, selfing = NULL, 
                            surv.matr = NULL, repr.matr = NULL, male.matr = NULL,
                            init.pop.sizes = NULL, num.gen = NULL, num.alleles = NULL, allele.freqs = NULL) {
 
 
rm(skeletonland)
skeletonland<-landscape.new.empty()
#define selfing rate
skeletonland<-landscape.new.floatparam(skeletonland,s=selfing)
#Hard coded in current generation, current epoch, max number generations, max number individuals
skeletonland<-landscape.new.intparam(skeletonland,h=num.pops,s=num.stgs,cg=0,ce=0,totgen=1000,maxland=20000)
#Hard coded in multiple paternity (yes) and density dependence (no) parameters
skeletonland<-landscape.new.switchparam(skeletonland,re=0,rd=0,mp=1,dd=0)

#local matrices, will give same demography to each local population
for (i in 1:num.pops)
	skeletonland<-landscape.new.local.demo(skeletonland,surv.matr, repr.matr, male.matr)

#cross habitat matrices
#epoch_s_matr<-matrix(0,nrow=4, ncol=4)
#epoch_r_matr<-matrix(0,nrow=4, ncol=4)
#epoch_m_matr<-matrix(0,nrow=4, ncol=4)

#no extinction allowed, hard coded
skeletonland<-landscape.new.epoch(skeletonland,epochprob=1,
    startgen=0,extinct=NULL,carry=carrying)

#LOCI
#Note that for SNPs, numalleles should be 2, allelesize only used for sequences
#type = 0 IAM, type = 1 SMM type = 2 DNA sequence
#assumes biparental transmission (transmission = 0)
rms.locus.type = NULL
if (loc.type == "SNP") {rms.locus.type = 2; num.alleles = 4; seq.length = rep(1,num.loc)}
if (loc.type == "microsat") rms.locus.type = 1
if (loc.type == "sequence") rms.locus.type = 2
for (l in 1:num.loc)
skeletonland<-landscape.new.locus(skeletonland, type=2, ploidy=2, mutationrate=mut.rate[l], 
    transmission=0, numalleles=num.alleles[l], frequencies=allele.freqs, allelesize=seq.length[l])

#assumes population initial sizes all defined nicely by user
skeletonland<-landscape.new.individuals(skeletonland,init.pop.sizes)

}
