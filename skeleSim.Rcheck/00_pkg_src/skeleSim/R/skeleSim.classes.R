#' @name skeleSim.classes
#' @import methods
setClassUnion("logOrNULL", c("logical", "NULL"))
setClassUnion("listOrNULL", c("list","NULL"))
setClassUnion("charOrNULL", c("character", "NULL"))
setClassUnion("intOrNum", c("integer","numeric", "NULL"))
setClassUnion("funcOrNULL", c("function", "NULL"))
setClassUnion("matrOrNULL", c("matrix", "NULL"))
setClassUnion("vectOrNULL", c("vector", "NULL"))
setClassUnion("posixOrNULL", c("POSIXct", "POSIXlt", "NULL"))

#' @rdname skeleSim.classes
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
#' @slot start.time a POSIXct representation of the starting time of the simulation.
#' @slot end.time a POSIXct representation of the end time of the simulation.
#' @slot num.reps number of replicates to run.
#' @slot sim.func a function that runs one replicate of the simulator.
#'   Must take and return only a \code{skeleSim.params} object.
#' @slot last.sample result of last call to \code{sim.func}.
#' @slot rep.analysis.func a function that analyzes the results of one
#'   simulation replicate.
#' @slot rep.result result from last call to \code{rep.analysis.func}.
#' @slot analysis.results a matrix containing result of all replicate analyses.
#' @slot sim.summary.func a function to summarize \code{rep.analysis}.
#' @slot summary.results a list containign result from call to
#'   \code{sim.summary.func}.
#' @slot sim.check.func a function to check the parameters object prior to
#'   running the simualtions
#' @slot sim.scen.checks a matrix containing results of 'checks' on scenario elements (T/F)
#' @slot other.checks a vector containing results of 'checks' on other param object elements
#' @slot scenario.reps a two column matrix describing which iteration matches
#'   which scenario/replicate
#' @slot analyses.requested vector of logicals specifying "Global", "Population",
#'   "Locus", or "Pairwise" analyses have been requested.
#'
setClass(
  Class = "skeleSim.params",
  slots = c(title = "charOrNULL", date = "posixOrNULL", quiet = "logOrNULL",
            question = "charOrNULL", simulator.type = "charOrNULL",
            simulator = "charOrNULL", wd = "charOrNULL",
            scenarios = "listOrNULL",
            num.reps = "intOrNum",  sim.func = "funcOrNULL",
            current.scenario = "intOrNum", current.replicate = "intOrNum",
            rep.sample = "ANY", rep.analysis.func = "funcOrNULL",
            rep.result = "intOrNum", analysis.results = "ANY",
            sim.summary.func = "funcOrNULL", summary.results = "listOrNULL",
            sim.check.func = "funcOrNULL", sim.scen.checks = "matrOrNULL",
            other.checks = "logOrNULL", scenario.reps = "matrOrNULL",
            analyses.requested = "logOrNULL"
  ),
  prototype = c(title = NULL, date = NULL, quiet = NULL, question = NULL,
                simulator.type = NULL, simulator = NULL, wd = NULL, scenarios = NULL,
                num.reps = NULL, sim.func = NULL,
                current.scenario = 1, current.replicate = NULL,
                rep.sample = NULL, rep.analysis.func = NULL, rep.result = NULL,
                analysis.results = NULL, sim.summary.func = NULL,
                summary.results = NULL, sim.check.func = NULL, sim.scen.checks = NULL,
                other.checks = NULL, scenario.reps = NULL,
                analyses.requested = c(Global = TRUE, Population = TRUE, Locus = TRUE,
                                       Pairwise = TRUE)
  )
)


#' @rdname skeleSim.classes
#' @title Scenario Parameters Class
#' @description An S4 class storing parameters for each simulation scenario
#'
#' @slot num.pops number of populations.
#' @slot pop.size a vector \code{num.pop} long giving size of each populaiton.
#' @slot sample.size a vector \code{num.pop} long giving the number of
#'   samples to take from each population.
#' @slot migration a \code{num.pop} x \code{num.pop} matrix giving the
#'   migration rates between each population.
#' @slot mig.helper a list of flags and values that are needed for the shiny interface but are not needed for the simulation
#'   itself.  Makes it easier to keep track of different ways to specify migration matrices for different scenarios.
#'   List elements will include migration model, rows and columns of landscape and distance function.
#' @slot simulator.params an object storing simulator-specific parameters. Can
#'   be a list or a simulator-specific class.
#'
setClass(
  Class = "scenario.params",
  slots = c(num.pops = "intOrNum", pop.size = "intOrNum",
            sample.size = "intOrNum", migration = "listOrNULL",
            mig.helper = "listOrNULL",
            locus.type = "charOrNULL", num.loci = "intOrNum",
            sequence.length = "intOrNum", mut.rate = "intOrNum",
            simulator.params = "ANY"
  ),
  prototype = c(num.pops = NULL, pop.size = NULL, sample.size = NULL,
      migration = NULL, mig.helper = NULL,
      locus.type = NULL, num.loci = NULL,
      sequence.length = NULL, mut.rate = NULL,
      simulator.params = NULL
      )
)
