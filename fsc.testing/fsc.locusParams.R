#' @name fscLocusParams
#' @title Create fastsimcoal locus parameter data.frames
#' @description Create fastsimcoal locus parameter data.frames
#'
#' @param sc a \linkS4class{scenario.params} object.
#' @param fsc a \linkS4class{fastsimcaoal.params} object.
#'
#' @return a \linkS4class{fastsimcoal.params} object with the
#'   \code{locus.params} slot filled.
#'
#' @export
#'
fsc.loadScenariolocusParams <- function(sc, fsc) {
  loci <- data.frame(
    type = sc@locus.type,
    num.loci = sc@num.loci,
    sequence.length = sc@sequence.length,
    mut.rate = sc@mut.rate
  )

  for(i in 1:nrow(loci)) {
    lp <- switch(loci$type[i],
      sequence = fsc.dna.locus(loci$sequence.length[i], loci$mut.rate[i]),
      microsat = fsc.msat.locus(loci$num.loci[i], loci$num.alleles[i], loci$mut.rate[i]),
      snp = fsc.snp.locus(loci$num.loci[i], loci$mut.rate[i])
    )
    loadLocusParams(fsc, lp)
  }
  fsc
}


#' @rdname fscLocusParams
#' @export
#'
fsc.locus.dna <- function(sequence.length, mut.rate, transition.rate = 1 / 3,
                          recomb.rate = 0, chromosome = 1) {
  df <- data.frame(
    chromosome = chromosome,
    type = "DNA",
    num.markers = sequence.length,
    recomb.rate = recomb.rate,
    param.4 = mut.rate,
    param.5 = transition.rate,
    param.6 = NA,
    stringsAsFactors = FALSE
  )
  df <- df[order(df$chromosome), ]
  attr(df, "ploidy") <- 1
  class(df) <- c(class(df), "fsc.locusParams")
  return(df)
}


#' @rdname fscLocusParams
#' @export
#'
fsc.locus.snp <- function(num.loci, min.freq, recomb.rate = 0, chromosome = 1) {
  df <- data.frame(
    chromosome = chromosome,
    type = "SNP",
    num.markers = num.loci,
    recomb.rate = recomb.rate,
    param.4 = min.freq,
    param.5 = NA,
    param.6 = NA,
    stringsAsFactors = FALSE
  )
  df <- df[order(df$chromosome), ]
  attr(df, "ploidy") <- 2
  class(df) <- c(class(df), "fsc.locusParams")
  return(df)
}


#' @rdname fscLocusParams
#' @export
#'
fsc.locus.msat <- function(num.loci, num.alleles, mut.rate, gsm.param = 0,
                           range.constraint = 0, recomb.rate = 0,
                           chromosome = 1) {
  df <- data.frame(
    chromosome = chromosome,
    type = "MICROSAT",
    num.markers = num.loci,
    recomb.rate = recomb.rate,
    param.4 = mut.rate,
    param.5 = gsm.param,
    param.6 = range.constraint,
    stringsAsFactors = FALSE
  )
  df <- df[order(df$chromosome), ]
  attr(df, "ploidy") <- 2
  class(df) <- c(class(df), "fsc.locusParams")
  return(df)
}
