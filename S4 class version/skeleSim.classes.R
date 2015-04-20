setClassUnion("logOrNULL", c("logical", "NULL"))
setClassUnion("listOrNULL", c("list","NULL"))
setClassUnion("charOrNULL", c("character", "NULL"))
setClassUnion("intOrNum", c("integer","numeric", "NULL"))
setClassUnion("intOrNULL", c("integer", "NULL"))
setClassUnion("funcOrNULL", c("function", "NULL"))
setClassUnion("posixOrNULL", c("POSIXct", "POSIXlt", "NULL"))


#' @title skeleSim Parameters Class
#' @description An S4 class storing generic parameters used throughout
#'   the workflow
#'
#' @rdname skeleSim.classes
#'
#' @slot title a title for the simulation. Used in labelling of output files.
#' @slot date datestamp for the simulation.
#' @slot quiet logical determining whether to limit progress reports.
#' @slot question a single character representing type of analytical question
#'   being addressed. Can be one of: (n)ull, p(o)wer, (p)erformance.
#' @slot simulator a single character representing which type of simulator
#'   to use. Can be one of: (c)oalescent or (f)orward-time.
#' @slot user.has.data logical recording whether or not the user has supplied
#'   data.
#' @slot scenarios a list of \code{scenario.params} objects.
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
#'
setClass(
  Class = "skeleSim.params",
  slots = c(title = "charOrNULL", date = "posixOrNULL", quiet = "logOrNULL",
            question = "charOrNULL", simulator = "charOrNULL",
            user.has.data = "logOrNULL", scenarios = "listOrNULL",
            num.reps = "intOrNULL", sim.func = "funcOrNULL",
            last.sample = "intOrNum", rep.analysis.func = "funcOrNULL",
            rep.result = "intOrNum", analysis.results = "intOrNum",
            sim.summary.func = "funcOrNULL", summary.results = "listOrNULL"
  ),
  prototype = c(title = NULL, date = NULL, quiet = NULL, question = NULL,
                simulator = NULL, user.has.data = NULL, scenarios = NULL,
                num.reps = NULL, sim.func = NULL, last.sample = NULL,
                rep.analysis.func = NULL, rep.result = NULL,
                analysis.results = NULL, sim.summary.func = NULL,
                summary.results = NULL
  )
)


#' @title Scenario Parameters Class
#' @description An S4 class storing parameters for each simulation scenario
#'
#' @rdname skeleSim.classes
#'
#' @slot num.pops number of populations.
#' @slot pop.size a vector \code{num.pop} long giving size of each populaiton.
#' @slot sample.size a vector \code{num.pop} long giving the number of
#'   samples to take from each population.
#' @slot migration a \code{num.pop} x \code{num.pop} matrix giving the
#'   migration rates between each population.
#' @slot simulator.params an object storing simulator-specific parameters. Can
#'   be a list or a simulator-specific class.
#'
setClass(
  Class = "scenario.params",
  slots = c(num.pops = "intOrNULL", pop.size = "intOrNULL",
            sample.size = "intOrNULL", migration = "intOrNum",
            simulator.params = "ANY"
  ),
  prototype = c(num.pops = NULL, pop.size = NULL, sample.size = NULL,
                migration = NULL, simulator.params = NULL
  )
)