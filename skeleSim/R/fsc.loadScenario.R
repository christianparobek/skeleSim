#' @title Load fastsimcoal locus parameters
#' @description Load fastsimcoal locus parameters
#'
#' @param num.pops number of populations.
#' @param pop.size a vector \code{num.pop} long giving size of each populaiton.
#' @param sample.size a vector \code{num.pop} long giving the number of
#'   samples to take from each population.
#' @param migration a \code{num.pop} x \code{num.pop} matrix giving the
#'   migration rates between each population.
#' @param locus.type a character representation of what type of marker to simulate.
#'   Can be "dna", "msat", or "snp".
#' @param num.loci number of msat or snp loci to simulate.
#' @param sequence.length number of DNA base pairs to use.
#' @param mut.rate mutation rate for DNA or msat.
#' @param sample.times a vector giving the number of generations in the past
#'   at which samples are taken.
#' @param growth.rate a vector giving the growth rate of each population.
#' @param hist.ev a matrix describing historical events.
#' @param num.chrom a value giving the number of chromosomes that the
#'   \code{locus.params} marker specifications should be copied for. If
#'   \code{NULL}, then chromosome assignment is taken from the
#'   \code{chromosome} column in \code{locus.params}. Any non-\code{NULL}
#'   integer will cause the \code{chromosome} column to be ignored.
#' @param transition.rate DNA: fraction of substitutions that are transitions.
#' @param recomb.rate recombination rate between adjacent markers.
#' @param chromosome number or character identifying which chromosome the marker
#'   is on.
#' @param min.freq SNP: minimum frequency for the derived allele.
#' @param gsm.param Value of the geometric parameter for a Generalized Stepwise
#'   Mutation (GSM) model. This value represents the proportion of mutations
#'   that will change the allele size by more than one step. Values between
#'   0 and 1 are required. A value of 0 is for a strict
#'   Stepwise Mutation Model (SMM).
#' @param range.constraint Range constraint (number of different alleles
#'   allowed). A value of 0 means no range constraint
#'
#' @export
#'
fsc.loadScenario <- function(num.pops, pop.size, sample.size, migration,
                             mut.rate, sample.times = NULL, growth.rate = NULL,
                             hist.ev = NULL, locus.type = c("dna", "msat", "snp"),
                             sequence.length = NULL, num.loci = NULL,
                             transition.rate = NULL, gsm.param = NULL,
                             range.constraint = NULL, min.freq = NULL,
                             recomb.rate = NULL, chromosome = NULL, num.chrom = NULL) {

  sc <- new("scenario.params")
  sc@num.pops <- num.pops
  sc@pop.size <- rep(pop.size, length.out = num.pops)
  sc@sample.size <- rep(sample.size, length.out = num.pops)
  sc@migration <- migration
  sc@locus.type <- locus.type
  sc@sequence.length <- sequence.length
  sc@num.loci <- num.loci
  sc@mut.rate <- mut.rate
  sc@migration <- migration

  fsc <- new("fastsimcoal.params")
  if(is.null(sample.times)) sample.times <- 0
  if(is.null(growth.rate)) growth.rate <- 0
  fsc@sample.times <- rep(sample.times, length.out = num.pops)
  fsc@growth.rate <- rep(growth.rate, length.out = num.pops)
  fsc@hist.ev <- hist.ev
  fsc@num.chrom <- num.chrom
  fsc@locus.params <- switch(match.arg(locus.type),
    dna = fsc.locus.dna(sequence.length, mut.rate, transition.rate, recomb.rate, chromosome),
    msat = fsc.locus.msat(num.loci, mut.rate, gsm.param, range.constraint, recomb.rate, chromosome),
    snp = fsc.locus.snp(num.loci, min.freq, recomb.rate, chromosome)
  )

  sc@simulator.params <- fsc
  return(sc)
}