setClassUnion("posixOrNULL", c("POSIXct", "POSIXlt", "NULL"))
setClassUnion("matrOrNULL", c("matrix", "NULL"))
setClassUnion("funcOrNULL", c("function", "NULL"))
setClassUnion("logOrNULL", c("logical", "NULL"))

#' @title skeleSim Parameters Class
#' @description An S4 class storing generic parameters used throughout
#'   the workflow
#'
#' @slot title a title for the simulation. Used in labelling of output files.
#' @slot date datestamp for the simulation.
#' @slot quiet logical determining whether to limit progress reports.
#' @slot question a single character representing type of analytical question
#'   being addressed. Can be one of: (n)ull, p(o)wer, (p)erformance.
#' @slot simulator.type a single character representing which type of simulator
#'   to use. Can be one of: (c)oalescent or (f)orward-time.
#' @slot simulator a three character code representing which simulator is being run.
#'   Currently codes for fastsimcoal(fsc) and rmetasim(rms) exist.
#' @slot scenarios a list of \code{scenario.params} objects.
#' @slot num.sim.reps number of replicates to run.
#' @slot sim.func a function that runs one replicate of the simulator.
#'   Must take and return only a \code{skeleSim.params} object.
#' @slot current.scenario number of current scenario being run.
#' @slot current.replicate number of current replicate within current scenario being run.
#' @slot rep.sample result of last call to \code{sim.func}.
#' @slot rep.analysis.func a function that analyzes the results of one
#'   simulation replicate.
#' @slot num.perm.reps number of permutation replicates to run for population structure
#'   statistics.
#' @slot rep.result result from last call to \code{rep.analysis.func}.
#' @slot analysis.results a matrix containing result of all replicate analyses.
#' @slot sim.summary.func a function to summarize \code{rep.analysis}.
#' @slot summary.results a list containign result from call to
#'   \code{sim.summary.func}.
#' @slot sim.check.func a function to check the parameters object prior to
#'   running the simualtions
#' @slot sim.scen.checks a matrix containing results of 'checks' on scenario elements (T/F)
#' @slot timing list containing elapsed time for a simulation
#' @slot other.checks a vector containing results of 'checks' on other param object elements
#' @slot scenario.reps a two column matrix describing which iteration matches
#'   which scenario/replicate
#' @slot analyses.requested vector of logicals specifying "Global", "Locus",
#'   or "Pairwise" analyses have been requested.
#'
#' @name skeleSim.classes
#' @aliases skeleSim.params skeleSim.params-class
#' @importFrom methods setClass new
#' @export
#'
skeleSim.params <- setClass(
  Class = "skeleSim.params",
  slots = c(
    title = "charOrNULL", date = "posixOrNULL", quiet = "logOrNULL",
    question = "charOrNULL", simulator.type = "charOrNULL",
    simulator = "charOrNULL", wd = "charOrNULL", scenarios = "listOrNULL",
    num.sim.reps = "intOrNum",  sim.func = "funcOrNULL",
    current.scenario = "intOrNum", current.replicate = "intOrNum",
    rep.sample = "ANY", rep.analysis.func = "funcOrNULL",
    num.perm.reps = "intOrNum", rep.result = "intOrNum",
    analysis.results = "ANY", sim.summary.func = "funcOrNULL",
    summary.results = "listOrNULL", sim.check.func = "funcOrNULL",
    sim.scen.checks = "matrOrNULL", other.checks = "logOrNULL",
    scenario.reps = "matrOrNULL", analyses.requested = "logOrNULL",
    timing = "listOrNULL"
  ),
  prototype = list(
    title = NULL, date = NULL, quiet = NULL, question = NULL,
    simulator.type = NULL, simulator = NULL, wd = NULL, scenarios = NULL,
    num.sim.reps = NULL, sim.func = NULL, current.scenario = 1,
    current.replicate = NULL, rep.sample = NULL, rep.analysis.func = NULL,
    num.perm.reps = NULL, rep.result = NULL,
    analysis.results = NULL, sim.summary.func = NULL, summary.results = NULL,
    sim.check.func = NULL, sim.scen.checks = NULL, other.checks = NULL,
    scenario.reps = NULL,
    analyses.requested = c(Global = TRUE, Locus = TRUE, Pairwise = TRUE),
    timing=NULL
  )
)


#' @slot num.pops number of populations.
#' @slot pop.size a vector \code{num.pop} long giving size of each population.
#' @slot sample.size a vector \code{num.pop} long giving the number of
#'   samples to take from each population.
#' @slot migration a list of one or more \code{num.pop} x \code{num.pop} matrices
#'   giving the migration rates between each population.
#' @slot locus.type a character representation of what type of marker to simulate.
#'   Can be "dna", "msat", or "snp".
#' @slot mig.helper a list of flags and values that are needed for the shiny interface but are not needed for the simulation
#'   itself.  Makes it easier to keep track of different ways to specify migration matrices for different scenarios.
#'   List elements will include migration model, rows and columns of landscape and distance function.
#' @slot num.loci number of msat or snp loci to simulate.
#' @slot sequence.length number of DNA base pairs to use.
#' @slot mut.rate mutation rate for DNA or msat.
#' @slot simulator.params an object storing simulator-specific parameters. Can
#'   be a list or a simulator-specific class.
#'
#' @rdname skeleSim.classes
#' @aliases scenario.params
#' @export
#'
scenario.params <- setClass(
  Class = "scenario.params",
  slots = c(
    num.pops = "intOrNum", pop.size = "intOrNum", sample.size = "intOrNum",
    migration = "listOrNULL", mig.helper = "listOrNULL",
    locus.type = "charOrNULL", num.loci = "intOrNum",
    sequence.length = "intOrNum", mut.rate = "intOrNum", simulator.params = "ANY"
  ),
  prototype = list(
    num.pops = NULL, pop.size = NULL, sample.size = NULL,
    migration = NULL,  mig.helper = NULL,
    locus.type = NULL, num.loci = NULL,
    sequence.length = NULL, mut.rate = NULL, simulator.params = NULL
  )
)
