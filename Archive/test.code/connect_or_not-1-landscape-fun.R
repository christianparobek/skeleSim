create.land <- function(
	num.pops= 2 ,
	num.stgs= 2 ,
	cg= 0 ,
	ce= 0 ,
	totgen= 1e+05 ,
	maxland= 1e+05 ,
	selfing= 0 ,
	re= 0 ,
	rd= 0 ,
	mp= 1 ,
	dd= 0 ,
	surv.matr = matrix(c(0.2, 0.4, 0, 0.1), nrow=2, ncol=2) ,
	repr.matr = matrix(c(0, 0, 4, 0), nrow=2, ncol=2) ,
	male.matr = matrix(c(0, 0, 0, 1), nrow=2, ncol=2) ,
	mig.rates = matrix(c(0, 0, 0, 0), nrow=2, ncol=2) ,
	carrying = c(2000, 100) ,
	mut.rate = c(1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04) ,
	num.alleles = c(5, 5, 5, 5, 5, 5, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8) ,
	init.pop.sizes = c(525, 525, 525, 525) ,
	allele.freqs = list(c(0.2, 0.2, 0.2, 0.2, 0.2), c(0.2, 0.2, 0.2, 0.2, 0.2), c(0.2, 0.2, 0.2, 0.2, 0.2), c(0.2, 0.2, 0.2, 0.2, 0.2), c(0.2, 0.2, 0.2, 0.2, 0.2), c(0.2, 0.2, 0.2, 0.2, 0.2), c(0.2, 0.2, 0.2, 0.2, 0.2), c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125), c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125), c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125), c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125), c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125), c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125), c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125), c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125), c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125), c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125), c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125), c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125), c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125))
)
{
skeletonland<-landscape.new.empty()
skeletonland<-landscape.new.floatparam(skeletonland,s=selfing)
skeletonland<-landscape.new.intparam(skeletonland,h=num.pops,s=num.stgs,cg=0,ce=0,totgen=100000,maxland=100000)
skeletonland<-landscape.new.switchparam(skeletonland,re=0,rd=0,mp=1,dd=0)
skeletonland<-landscape.new.local.demo(skeletonland,surv.matr, repr.matr, male.matr)
epoch_s_matr<-matrix(0,nrow=(num.pops*num.stgs), ncol=(num.pops*num.stgs))
epoch_r_matr <- landscape.mig.matrix(h=num.pops,s=num.stgs,R.custom=mig.rates)$R
epoch_m_matr <- epoch_s_matr
skeletonland<-landscape.new.epoch(skeletonland,epochprob=1, epoch_s_matr, epoch_r_matr, epoch_m_matr,startgen=0,extinct=NULL,carry=carrying)
skeletonland<-landscape.new.locus(skeletonland, type=1, ploidy=2, mutationrate=mut.rate[1], transmission=0, numalleles=num.alleles[1], frequencies=allele.freqs[[1]])
skeletonland<-landscape.new.locus(skeletonland, type=1, ploidy=2, mutationrate=mut.rate[2], transmission=0, numalleles=num.alleles[2], frequencies=allele.freqs[[2]])
skeletonland<-landscape.new.locus(skeletonland, type=1, ploidy=2, mutationrate=mut.rate[3], transmission=0, numalleles=num.alleles[3], frequencies=allele.freqs[[3]])
skeletonland<-landscape.new.locus(skeletonland, type=1, ploidy=2, mutationrate=mut.rate[4], transmission=0, numalleles=num.alleles[4], frequencies=allele.freqs[[4]])
skeletonland<-landscape.new.locus(skeletonland, type=1, ploidy=2, mutationrate=mut.rate[5], transmission=0, numalleles=num.alleles[5], frequencies=allele.freqs[[5]])
skeletonland<-landscape.new.locus(skeletonland, type=1, ploidy=2, mutationrate=mut.rate[6], transmission=0, numalleles=num.alleles[6], frequencies=allele.freqs[[6]])
skeletonland<-landscape.new.locus(skeletonland, type=1, ploidy=2, mutationrate=mut.rate[7], transmission=0, numalleles=num.alleles[7], frequencies=allele.freqs[[7]])
skeletonland<-landscape.new.locus(skeletonland, type=1, ploidy=2, mutationrate=mut.rate[8], transmission=0, numalleles=num.alleles[8], frequencies=allele.freqs[[8]])
skeletonland<-landscape.new.locus(skeletonland, type=1, ploidy=2, mutationrate=mut.rate[9], transmission=0, numalleles=num.alleles[9], frequencies=allele.freqs[[9]])
skeletonland<-landscape.new.locus(skeletonland, type=1, ploidy=2, mutationrate=mut.rate[10], transmission=0, numalleles=num.alleles[10], frequencies=allele.freqs[[10]])
skeletonland<-landscape.new.locus(skeletonland, type=1, ploidy=2, mutationrate=mut.rate[11], transmission=0, numalleles=num.alleles[11], frequencies=allele.freqs[[11]])
skeletonland<-landscape.new.locus(skeletonland, type=1, ploidy=2, mutationrate=mut.rate[12], transmission=0, numalleles=num.alleles[12], frequencies=allele.freqs[[12]])
skeletonland<-landscape.new.locus(skeletonland, type=1, ploidy=2, mutationrate=mut.rate[13], transmission=0, numalleles=num.alleles[13], frequencies=allele.freqs[[13]])
skeletonland<-landscape.new.locus(skeletonland, type=1, ploidy=2, mutationrate=mut.rate[14], transmission=0, numalleles=num.alleles[14], frequencies=allele.freqs[[14]])
skeletonland<-landscape.new.locus(skeletonland, type=1, ploidy=2, mutationrate=mut.rate[15], transmission=0, numalleles=num.alleles[15], frequencies=allele.freqs[[15]])
skeletonland<-landscape.new.locus(skeletonland, type=1, ploidy=2, mutationrate=mut.rate[16], transmission=0, numalleles=num.alleles[16], frequencies=allele.freqs[[16]])
skeletonland<-landscape.new.locus(skeletonland, type=1, ploidy=2, mutationrate=mut.rate[17], transmission=0, numalleles=num.alleles[17], frequencies=allele.freqs[[17]])
skeletonland<-landscape.new.locus(skeletonland, type=1, ploidy=2, mutationrate=mut.rate[18], transmission=0, numalleles=num.alleles[18], frequencies=allele.freqs[[18]])
skeletonland<-landscape.new.locus(skeletonland, type=1, ploidy=2, mutationrate=mut.rate[19], transmission=0, numalleles=num.alleles[19], frequencies=allele.freqs[[19]])
skeletonland<-landscape.new.locus(skeletonland, type=1, ploidy=2, mutationrate=mut.rate[20], transmission=0, numalleles=num.alleles[20], frequencies=allele.freqs[[20]])
skeletonland<-landscape.new.individuals(skeletonland,init.pop.sizes)
}