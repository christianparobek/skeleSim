rm(list = ls())
source("skeleSim.classes.R")
source("fastsimcoal.skeleSim.R")
source("skeleSim.funcs.R")

# ---- Load parameters ---

# create new skeleSim parameters object
test.params <- new("skeleSim.params")
test.params@title <- "testRun"
test.params@date <- Sys.time()
test.params@quiet <- FALSE
test.params@question <- "n"
test.params@simulator.type <- "c"
test.params@simulator <- "fsc"
test.params@num.reps <- 100
test.params@timing <- 2
test.params@sim.func <- fsc.run
test.params@wd <- "testRun.wd"

# create a base scenario parameters object. It will be copied and modified
#   later for different scenarios
base.scenario <- new("scenario.params")
base.scenario@num.pops <- 3
base.scenario@pop.size <- c(50, 100, 500)
base.scenario@sample.size <- c(25, 50, 25)
base.scenario@migration <- list(matrix(
  c(0, 0.01, 0.05, 0.025, 0, 0.025, 0.05, 0.01, 0),
  nrow = base.scenario@num.pops
))
base.scenario@locus.type <- "sequence"
base.scenario@num.loci <- 1
base.scenario@sequence.length <- 400
base.scenario@mut.rate <- 1e-4

# create fastsimcoal params object to load into base scenario
fsc.params <- new("fastsimcoal.params")
# to change the executable, either explicitly set the @fastsimcoal.exec slot or
# initialize with:
#   fsc.params <- new("fastsimcoal.params", fastsimcoal.exec = "fsc252")
fsc.params@sample.times <- c(0, 0, 0)
fsc.params@growth.rate <- c(0, 0, 0)
fsc.params@hist.ev <- fsc.histEvMat(0)
#fsc.params@hist.ev[, "num.gen"] <- 100
#fsc.params@hist.ev[, "sink.deme"] <- 1
fsc.params@locus.params <- fsc.locusParamsMat(base.scenario)
fsc.params@inf.site.model <- FALSE
# load fastsimcoal params
base.scenario@simulator.params <- fsc.params

# create list of scenarios and modify
scenario.list <- lapply(1:3, function(i) base.scenario)
#  decrease the mutation rate in scenario 2...
scenario.list[[2]]@mut.rate <- 1e-5
scenario.list[[2]]@simulator.params@locus.params <- fsc.locusParamsMat(scenario.list[[2]])
#  decrease the migration rate in scenario 3...
scenario.list[[3]]@migration[[1]] <- scenario.list[[3]]@migration[[1]] * 0.1

# load scenarios
test.params@scenarios <- scenario.list


# ---- Set analysis function ----
test.params@rep.analysis.func <- function(params) {
  result = rnorm(5)
  names(result) <- paste("result", 1:5, sep = ".")
  params@rep.result <- result
  params
}

# --- Set parameter check function ---
test.params@sim.check.func <- fsc.scenarioCheck


# ---- Run replicates ----
test.params <- runSim(test.params)


# ---- Summarize analysis results ----
test.params <- summ.stats.table(test.params)
