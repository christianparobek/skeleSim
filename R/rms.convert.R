#' @title Extract genetic data from landscape object
#' @description Extract genetic data from landscape object
#'
#' @param Rland a \code{rmetasim} landscape object.
#' @param locus.type character giving type of locus (\code{microsat} or
#'   \code{sequence})
#'
#' @return A gtypes object containing genotypes (or sequences)
#'
#' @importFrom rmetasim landscape.populations is.landscape landscape.make.genind landscape.locus.states
#' @importFrom ape as.DNAbin
#' @importFrom strataG genind2gtypes sequence2gtypes
#'
#' @export
#'
rms.convert <- function(Rland, locus.type) {
  ltype <- locus.type

  # print("in rms.convert")
  # print(paste("locustype:",ltype))

  if (!is.landscape(Rland)) stop("incoming landscape problem in rms.convert")

  this.rep.result=NULL
  if (ltype%in%c("microsatellite","MICROSAT","microsat")) {
    this.rep.result <- genind2gtypes(landscape.make.genind(Rland))
  } else if (ltype=="sequence") {
    # print("converting rmetasim sequences")
    states <- as.data.frame(landscape.locus.states(Rland,1))
    genos <- data.frame(pop=landscape.populations(Rland),aindex=Rland$individuals[,7])
    seq <- merge(genos,states,all.x=T)
    seq <- seq[order(seq$pop),]
    dna.seq <- strsplit(as.character(tolower(seq$state)),"")
    dnabin <- as.DNAbin(do.call(rbind,strsplit(tolower(as.character(seq$state)),"")))
    this.rep.result <- sequence2gtypes(strata=seq$pop,x=dnabin)
  } else if (ltype=="SNP") {
    this.rep.result <- genind2gtypes(landscape.make.genind(Rland))
    #this.rep.result <- genind2gtypes(landscape.make.genind(landscape.snp.convert(Rland)))
  }

  this.rep.result
}
