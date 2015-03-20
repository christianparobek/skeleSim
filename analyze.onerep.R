##### Function to take mega list,
##    1) detect which type: genind, or DNAbin plus strata file, if DNAbin: convert to strataG (gtypes) with DNAbin2gtypes()
##    2) check for population structure, Y/N
##    3) run all for that type of data  #DNAbin or genind

function(mega.list, FUN){
  if(class(mega.list$rep.result) == "genind"){
    if(strata < 1){
      ## haplotype diversity and % unique haps per population
      ## strata.split is a bit slow for more populations
      df.hap.div  <- summary(mega.list$rep.result)$by.strata[,3]
      df.pct.haps <- summary(mega.list$rep.result)$by.strata[,4]

      if(gtypeobj )


    } else if (type = adegenit-SNPs, msat) {
      if(ploidy(mega.list$rep.result) ==
           pop <- strata
    }

  }

}
