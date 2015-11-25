#' @title Extract genetic data from landscape object
#' @description Extract genetic data from landscape object
#'
#' @param Rland a \code{rmetasim} landscape object.
#' @param locus.type character giving type of locus (\code{microsat} or
#'   \code{sequence})
#'
#' @importFrom rmetasim landscape.populations
#'
rms.convert <- function(Rland, locus.type) {
        ltype <- locus.type
        if (ltype=="microsat")
            {
                this.rep.result <- landscape.make.genind(Rland)
            }
        else if (ltype=="sequence")
            {
                states <- as.data.frame(landscape.locus.states(1,Rland))
                genos <- data.frame(pop=landscape.populations(Rland),aindex=Rland$individuals[,7])
                seq <- merge(genos,states,all.x=T)
                seq <- seq[order(seq$pop),]
                dna.seq <- strsplit(as.character(tolower(seq$state)),"")
                this.rep.result <- list(strata=data.frame(seq$pop),
                                           dna.seq=as.DNAbin(do.call(rbind,strsplit(tolower(as.character(seq$state)),""))))
            }
        else if (ltype=="SNP")
            {
                this.rep.result <- landscape.make.genind(landscape.snp.convert(Rland))
            }
      this.rep.result

}
