#' @title Create genind object from landscape
#' @description Create genind object from landscape
#'
#' @param l a \code{rmetasim} landscape object.
#' @param popnames character vector of population names.
#'
#' @importFrom adegenet genind
#'
#' @export
#' 
landscape.make.genind <- function(Rland,popnames=NULL)
  {
      tab <- 2*landscape.ind.freq(Rland)
      dimnames(tab) <- list(rownames=1:dim(tab)[1],colnames=landscape.freq.locnames(Rland))
      if (is.null(popnames))
          {
              populations <- landscape.populations(Rland)
          }
      else
          {
              populations <- popnames
          }
      adegenet::genind(tab,pop=as.factor(populations),ploidy=2)
  }
