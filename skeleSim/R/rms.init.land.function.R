#' @title Initialize a landscape object
#' @description Initialize a landscape object
#'
#' @param num.pops ?
#' @param carrying ?
#' @param sample.size ?
#' @param mig.rates ?
#' @param num.loc ?
#' @param loc.type ?
#' @param mut.rate ?
#' @param seq.length ?
#' @param num.stgs ?
#' @param selfing ?
#' @param surv.matr ?
#' @param repr.matr ?
#' @param male.matr ?
#' @param init.pop.sizes ?
#' @param num.gen ?
#' @param num.alleles ?
#' @param allele.freqs ?
#'
#' @importFrom rmetasim landscape.new.empty landscape.new.floatparam
#'   landscape.new.intparam landscape.new.switchparam landscape.new.local.demo
#'   landscape.mig.matrix landscape.new.epoch landscape.new.locus landscape.new.individuals
#' @export
rms.init.landscape.func <- function(num.pops = NULL, carrying = NULL, sample.size = NULL,
                                    mig.rates = NULL, num.loc = NULL, loc.type = NULL,
                                    mut.rate = NULL, seq.length = NULL, num.stgs = NULL,
                                    selfing = NULL, surv.matr = NULL, repr.matr = NULL,
                                    male.matr = NULL, init.pop.sizes = NULL,
                                    num.gen = NULL, num.alleles = NULL, allele.freqs = NULL)
{

matstrfun <- function(mat)
{
    if (!is.null(mat))
        {
            ms <- "matrix(c("
            ms <- c(ms,paste(mat,collapse=", "))
            ms <- c(ms,paste0("), nrow=",dim(mat)[1],", ncol=",dim(mat)[2],")"))
            paste(ms,collapse="")
        } else {NULL}
}
vecstrfun <- function(vec)
{
    if (!is.null(vec))
        {
            ms <- "c("
            ms <- c(ms,paste(vec,collapse=", "))
            ms <- c(ms,")")
            paste(ms,collapse="")
        } else {NULL}
}
listvecstrfun <- function(l)
    {
        if (!is.null(l))
            {
                paste0("list(",paste0(sapply(l,vecstrfun),collapse=", "),")")
            } else {NULL}
    }

    fstr <- c("create.land <- function(",
              paste("num.pops=",num.pops,","),
              paste("num.stgs=",num.stgs,","),
              paste("cg=",0,","),
              paste("ce=",0,","),
              paste("totgen=",100000,","),
              paste("maxland=",100000,","),
              paste("selfing=",selfing,","),
              paste("re=",0,","),
              paste("rd=",0,","),
              paste("mp=",1,","),
              paste("dd=",0,","),
              paste("surv.matr =",matstrfun(surv.matr),","),
              paste("repr.matr =",matstrfun(repr.matr),","),
              paste("repr.matr =",matstrfun(male.matr),","),
              paste("mig.rates =",matstrfun(mig.rates),","),
              paste("carrying =",vecstrfun(carrying),","),
              paste("mut.rate =",vecstrfun(mut.rate),","),
              paste("num.alleles =",vecstrfun(num.alleles),","),
              paste("allele.freqs =",listvecstrfun(allele.freqs))
              )
    fstr <- c(fstr,")\n{")
              
fstr <- c(fstr,"skeletonland<-landscape.new.empty()")
#define selfing rate
fstr <- c(fstr,
          paste0("skeletonland<-landscape.new.floatparam(skeletonland,s=selfing)"))
###Hard coded in current generation, current epoch, max number generations, max number individuals
fstr <- c(fstr,paste0("skeletonland<-landscape.new.intparam(skeletonland,h=num.pops,s=num.stgs,cg=0,ce=0,totgen=100000,maxland=100000)"))

###Hard coded in multiple paternity (yes) and density dependence (no) parameters
fstr <- c(fstr,paste0("skeletonland<-landscape.new.switchparam(skeletonland,re=0,rd=0,mp=1,dd=0)"))

#local matrices, will give same demography to each local population
fstr <- c(fstr,paste0("skeletonland<-landscape.new.local.demo(skeletonland,surv.matr, repr.matr, male.matr)"))

#cross habitat matrices

fstr <- c(fstr,paste0("epoch_s_matr<-matrix(0,nrow=(num.pops*num.stgs), ncol=(num.pops*num.stgs))"))
fstr <- c(fstr,paste0("epoch_r_matr <- landscape.mig.matrix(h=num.pops,s=num.stgs,R.custom=mig.rates)$R"))
fstr <- c(fstr,paste0("epoch_m_matr <- epoch_s_matr"))
    
#no extinction allowed, hard coded
fstr <- c(fstr,paste0("skeletonland<-landscape.new.epoch(skeletonland,epochprob=1, epoch_s_matr, epoch_r_matr, epoch_m_matr,startgen=0,extinct=NULL,carry=carrying)"))

#LOCI
#Note that for SNPs, numalleles should be 2, allelesize only used for sequences
#type = 0 IAM, type = 1 SMM type = 2 DNA sequence
#assumes biparental transmission (transmission = 0)
rms.locus.type = NULL

    print(loc.type)
    
if (loc.type == "SNP") {rms.locus.type = 2; num.alleles = 4; seq.length = rep(1,num.loc)}
if (loc.type %in% c("microsat","MICROSAT","microsatellite")) rms.locus.type = 1
if (loc.type == "sequence") rms.locus.type = 2
    for (l in 1:num.loc)
        if (rms.locus.type==2)
        {
            if (l==1) #only one sequence locus possible and it creates a maternally inherited haploid locus
                fstr <- c(fstr,paste0("skeletonland<-landscape.new.locus(skeletonland, type=2, ploidy=1, mutationrate=mut.rate[",l,"]",
                                      "transmission=1, numalleles=num.alleles[",l,"], frequencies=,allele.freqs[[",l,"]])")) #temp remove allele size
            
        }
        else
        {
            fstr <- c(fstr,paste0("skeletonland<-landscape.new.locus(skeletonland, type=1, ploidy=2, mutationrate=mut.rate[",l,"], transmission=0, numalleles=num.alleles[",l,"], frequencies=allele.freqs[[",l,"]])"))
        }
#assumes population initial sizes all defined nicely by user
    fstr <- c(fstr,"skeletonland<-landscape.new.individuals(skeletonland,init.pop.sizes)")
    fstr <- c(fstr,"}")
paste(fstr,collapse="\n")
}

