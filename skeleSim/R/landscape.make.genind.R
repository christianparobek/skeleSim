landscape.make.genind <-
function(l,popnames=NULL)
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
