#' @title Extract genetic data from landscape object
#' @description Extract genetic data from landscape object
#'
#' @param Rland a \code{rmetasim} landscape object.
#' @param locus.type character giving type of locus (\code{microsat} or
#'   \code{sequence})
#'
#' @importFrom rmetasim landscape.populations is.landscape landscape.make.genind
#' @export
rms.convert <- function(Rland, locus.type) {
    
    ltype <- locus.type

    print("in rms.convert")
    print(paste("locustype:",ltype))
    
    if (!is.landscape(Rland)) {stop("incoming landscape problem in rms.convert")}
    
    
    
    this.rep.result=NULL
    if (ltype%in%c("microsatellite","MICROSAT","microsat"))
    {
        this.rep.result <- landscape.make.genind(Rland)
    }
    else if (ltype=="sequence")
    {
        print("converting rmetasim sequences")
        states <- as.data.frame(landscape.locus.states(Rland,1))
        genos <- data.frame(pop=landscape.populations(Rland),aindex=Rland$individuals[,7])
        seq <- merge(genos,states,all.x=T)
        seq <- seq[order(seq$pop),]
        dna.seq <- strsplit(as.character(tolower(seq$state)),"")
        this.rep.result <- list(strata=data.frame(seq$pop),
                                dna.seq=new("multidna",as.DNAbin(do.call(rbind,strsplit(tolower(as.character(seq$state)),"")))))
    }
    else if (ltype=="SNP")
    {
        this.rep.result <- landscape.make.genind(landscape.snp.convert(Rland))
    }
    

    this.rep.result
}
