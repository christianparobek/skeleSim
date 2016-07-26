#' @title  Garza Williamson M Ratio (bottleneck) Statistic, Per Pop, Per Loc
#' @description calculates the mratio from a gtypes object
#' @param gen.data.gtype a gtype object
#' @export
calc.mratio <- function(gen.data.gtype) {
  pop.locus.df <- as.matrix(expand.grid(pop = 1:length(unique(gen.data.gtype@strata)),
                                        locus = 1:length(gen.data.gtype@loci)))

  calc.mratio.2<- function(all.freq){
    #find the smallest allele present
    if (any(all.freq >0)){
      min.all <- min(which(all.freq>0))
      #find the largest allele present
      max.all <- max(which(all.freq>0))
      #calculate the number of alleles present
      sum.all.present <- sum(all.freq>0)
      sum.all.present/(max.all-min.all+1)
    }
    else return(NA)
  }

  all_freq_list<-alleleFreqs(gen.data.gtype, by.strata=T)
  mrat_pop_and_loc<-apply(pop.locus.df,1, function(x) {
    calc.mratio.2(all_freq_list[[x[2]]][,,x[1]][,1])
  })

  all_freq_list<-alleleFreqs(gen.data.gtype, by.strata=F)
  mrat_loc<-sapply(1:length(all_freq_list), function(x){
    calc.mratio.2(as.vector(all_freq_list[[x]][,1]))
  })
  return(c(mrat_loc,mrat_pop_and_loc))

}

