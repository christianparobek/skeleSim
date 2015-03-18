
#
#
#
# function to convert genotypes into frequencies per individual
#
landscape.ind.freq <- function(l)
  {
    aml <- vector("list",length(landscape.ploidy(l)))
    for (loc in 1:length(landscape.ploidy(l)))
      {
        genos <- landscape.locus(loc,l)[,-1:-landscape.democol()]
        ploidy <- landscape.ploidy(l)[loc]
        amat <- sapply(as.numeric(names(table(genos))),function(x,genos,pl)
                       {
                         if (pl==2)
                           {
                             (genos[,1]==x)+(genos[,2]==x)
                           } else
                         {
                           genos==x
                         }
                       },genos=genos,pl=ploidy)
        aml[[loc]] <- apply(amat,2,function(x,pl){x/pl},pl=ploidy) #allele freqs per ind
      }
    do.call(cbind,aml)
  }

#
#
#
# function to provide locus/allele names
#
#
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



landscape.make.genind <- function(l,popnames=NULL)
  {
      require(adegenet)
      
      tab <- 2*landscape.ind.freq(l)
      dimnames(tab) <- list(rownames=1:dim(tab)[1],colnames=landscape.freq.locnames(l))
      if (is.null(popnames))
          {
              populations <- landscape.populations(l)
          }
      else
          {
              populations <- popnames
          }
      genind(tab,pop=as.factor(populations),ploidy=2)
  }

landscape.make.genpop <- function(l)
    {
        genind2genpop(landscape.make.genind(l))
    }
