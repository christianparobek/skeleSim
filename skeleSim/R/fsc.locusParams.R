#' @name fsc.locusParams
#' @title Create fastsimcoal locus parameter data.frames
#' @description Create fastsimcoal locus parameter data.frames
#'
#' @param sequence.length DNA: number of base pairs for sequence.
#' @param mut.rate DNA: mutation rate per bp, MICROSAT: mutation rate per locus.
#' @param transition.rate DNA: fraction of substitutions that are transitions.
#' @param recomb.rate recombination rate between adjacent markers.
#' @param chromosome number or character identifying which chromosome the marker
#'   is on.
#' @param num.loci MICROSAT, SNP: number of loci to simulate.
#' @param min.freq SNP: minimum frequency for the derived allele.
#' @param gsm.param Value of the geometric parameter for a Generalized Stepwise
#'   Mutation (GSM) model. This value represents the proportion of mutations
#'   that will change the allele size by more than one step. Values between
#'   0 and 1 are required. A value of 0 is for a strict
#'   Stepwise Mutation Model (SMM).
#' @param range.constraint Range constraint (number of different alleles
#'   allowed). A value of 0 means no range constraint
#'
NULL

#' @rdname fsc.locusParams
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


#' @rdname fsc.locusParams
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


#' @rdname fsc.locusParams
#' @export
#'
fsc.locus.msat <- function(num.loci, mut.rate, gsm.param = 0,
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