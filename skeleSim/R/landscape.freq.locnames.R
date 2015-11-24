landscape.freq.locnames <- function(l)
  {
    namevec <- NULL
    for (loc in 1:length(landscape.ploidy(l)))
      {
        genos <- landscape.locus(loc,l)[,-1:-landscape.democol()]
        namevec <- c(namevec,paste("L",paste(loc,names(table(unlist(genos))),sep="."),sep=''))
      }
    namevec
  }
