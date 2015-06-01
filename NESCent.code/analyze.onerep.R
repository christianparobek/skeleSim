##### Function to take mega list,
##    1) detect which type: genind, or DNAbin plus strata file, if DNAbin: convert to strataG (gtypes) with DNAbin2gtypes()
##    2) check for population structure, Y/N
##    3) run all for that type of data  #DNAbin or genind

analyze.onerep <- function(mega.list, FUN){
  if(class(mega.list$rep.result) == "genind"){

    ## split by strata needs to be
    xhierfstat <- genind2hierfstat(mega.list$rep.result, pop=mega.list$rep.result$pop) # @ or $ pop?  mega.list$rep.result@pop
    Summary<-basic.stats(xhierfstat,diploid=TRUE,digits=4)#needs to be separated, so that it gets analysed per population
    Hs<-Summary$overall[c(2)]#this pulls out gene diversity, Hs
    Ho<-Summary$overall[c(1)]#this pulls out obeserved heterozygosities, Ho
    Fis<-Summary$overall[c(9)]#this pulls out Fis, following Nei(1987)

  } else {

    gtype.rep.result <- DNAbin2gtypes(mega.list$rep.result, mega.list$rep.result$strata, ) # From function DNAbin2gtypes.R, read.FASTA(), and strata

    if(gtype.rep.result strata > 1){ ## ONLY NEEDED FOR SLOW OF strata.split
      ## haplotype diversity and % unique haps per population
      ## strata.split is a bit slow for more populations
      df.hap.div  <- summary(mega.list$rep.result)$by.strata[,3]
      df.pct.haps <- summary(mega.list$rep.result)$by.strata[,4]

      if(gtypeobj )
    } else {
      ## ARe there other non-split needing by population ones... maybe get rid of that if else

    }

    } else if (type = adegenit-SNPs, msat) {
      if(ploidy(mega.list$rep.result) ==
           pop <- strata
    }

  }

#THIS IS TO NAME VECTOR ELEMENTS FOR THE SUMMARY STATS
names(summarystatvectorNAME) = paste("summarystatvectorNAME",1:length(summarystatvector),sep="_")

#TO CONCATENATE THE NAMED SUMMARY STAT VECTORS, to be passed to Christian
c(summarystatvectorNAME,summarystatvectorNAME1,....)


}
