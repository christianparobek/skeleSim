#' @title Load skeleSim scenario parameters for fastsimcoal
#' @description Load skeleSim scenario parameters for fastsimcoal
#'
#' @param num.pops number of populations.
#' @param pop.size a vector giving size of each populaiton.
#' @param sample.size a vector giving the number of samples to take from each
#'   population.
#' @param migration a \code{num.pop} x \code{num.pop} matrix or list of matrices
#'   giving the migration rates between each population.
#' @param locus.type a character representation of what type of marker to simulate.
#'   Can be "dna", "msat", or "snp".
#' @param num.loci \code{msat, snp}: number of loci to simulate.
#' @param sequence.length \code{dna}: number of DNA base pairs to use.
#' @param mut.rate \code{dna, msat}: per base pair or locus mutation rate.
#' @param sample.times a vector giving the number of generations in the past
#'   at which samples are taken.
#' @param growth.rate a vector giving the growth rate of each population.
#' @param hist.ev a matrix describing historical events.
#' @param num.chrom a value giving the number of chromosomes that the
#'   \code{locus.params} marker specifications should be copied for. If
#'   \code{NULL}, then chromosome assignment is taken from the
#'   \code{chromosome} column in \code{locus.params}. Any non-\code{NULL}
#'   integer will cause the \code{chromosome} column to be ignored.
#' @param transition.rate dna: fraction of substitutions that are transitions.
#' @param recomb.rate recombination rate between adjacent markers.
#' @param chromosome number or character identifying which chromosome the marker
#'   is on.
#' @param min.freq \code{snp}: minimum frequency for the derived allele.
#' @param gsm.param \code{msat}: Value of the geometric parameter for a
#'   Generalized Stepwise Mutation (GSM) model. This value represents the
#'   proportion of mutations that will change the allele size by more than
#'   one step. Values between 0 and 1 are required. A value of 0 is for a
#'   strict Stepwise Mutation Model (SMM).
#' @param range.constraint \code{msat}: Range constraint (number of different
#'   alleles allowed). A value of 0 means no range constraint.
#'
#' @note Vectors for \code{pop.size, sample.size, sample.times, and growth.rate}
#'   will be expanded/recycled to ensure they are as long as \code{num.pops}.
#'
#' Depending on the choice of \code{locus.type}, values for some arguments may
#'   be ignored. See argument list above for which arguments are applicable
#'   to which \code{locus.type}.
#'
#'
#'
#' @return a \linkS4class{scenario.params} object to be loaded into a list in the
#'   \code{scenarios} slot of a \linkS4class{skeleSim.params} object.
#'
#' @export
#'
fsc.loadScenario <- function(
  num.pops, pop.size, sample.size, mut.rate, migration = NULL, sample.times = NULL,
  growth.rate = NULL, hist.ev = NULL, locus.type = c("dna", "msat", "snp"),
  sequence.length = NULL, num.loci = NULL, transition.rate = NULL, gsm.param = NULL,
  range.constraint = NULL, min.freq = NULL, recomb.rate = NULL, chromosome = NULL,
  num.chrom = NULL) {

  # convert NULLs to default values (some may not be used)
  if(is.null(num.loci)) num.loci <- 1
  if(is.null(sample.times)) sample.times <- 0
  if(is.null(growth.rate)) growth.rate <- 0
  if(is.null(transition.rate)) transition.rate <- 1 / 3
  if(is.null(gsm.param)) gsm.param <- 0
  if(is.null(range.constraint)) range.constraint <- 0
  if(is.null(recomb.rate)) recomb.rate <- 0
  if(is.null(chromosome)) chromosome <- 1

  pop.size <- rep(pop.size, length.out = num.pops)
  sample.size <- rep(sample.size, length.out = num.pops)

  # load general scenario parameters
  sc <- new("scenario.params")
  sc@num.pops <- num.pops
  sc@pop.size <- pop.size
  sc@sample.size <- sample.size
  sc@locus.type <- locus.type
  sc@sequence.length <- sequence.length
  sc@num.loci <- num.loci
  sc@mut.rate <- mut.rate
  sc@migration <- if(is.matrix(migration)) {
    list(migration)
  } else if(is.list(migration)) {
    migration
  } else NULL

  # load fastsimcoal-specific parameters
  fsc <- new("fastsimcoal.params")
  fsc@pop.info = fscPopInfo(
    pop.size = pop.size,
    sample.size = sample.size,
    sample.times = rep(sample.times, length.out = num.pops),
    growth.rate = rep(growth.rate, length.out = num.pops)
  )
  fsc@hist.ev <- hist.ev
  fsc@locus.params <- fscLocusParams(
    locus.type = locus.type,
    sequence.length = sequence.length,
    num.loci = num.loci,
    mut.rate = mut.rate,
    transition.rate = transition.rate,
    gsm.param = gsm.param,
    range.constraint = range.constraint,
    recomb.rate = recomb.rate,
    chromosome = chromosome,
    num.chrom = num.chrom,
    ploidy = ploidy
  )

  sc@simulator.params <- fsc
  return(sc)
}