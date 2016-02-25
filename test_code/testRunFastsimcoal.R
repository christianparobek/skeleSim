rm(list = ls())
library(skeleSim)

# ---- Load parameters ---

# create new skeleSim parameters object
test.params <- new("skeleSim.params")
test.params@title <- "testRun"
test.params@date <- Sys.time()
test.params@quiet <- FALSE
test.params@question <- "n"
test.params@simulator.type <- "c"
test.params@simulator <- "fsc"
test.params@num.reps <- 10
test.params@num.perm.reps <- 100
test.params@num.cores <- 1
test.params@sim.func <- fsc.run
test.params@wd <- "testRun.wd"

# create a base scenario parameters object. It will be copied and modified
#   later for different scenarios
base.scenario <- fsc.loadScenario(
  num.pops = 3,
  pop.size = c(50, 100, 500),
  sample.size = c(25, 50, 10),
  migration = list(matrix(
    c(0, 0.01, 0.05, 0.025, 0, 0.025, 0.05, 0.01, 0), nrow = 3
  )),
  locus.type = "msat",
  num.loci = 20,
  mut.rate = 1e-2,
  range.constraint = c(0, 20)
)

# create list of scenarios and modify
scenario.list <- lapply(1:3, function(i) base.scenario)
#  decrease the mutation rate in scenario 2...
scenario.list[[2]]@mut.rate <- 1e-5
locus.params <- fsc.locus.dna(base.scenario@sequence.length, scenario.list[[2]]@mut.rate)
scenario.list[[2]]@simulator.params@locus.params <- locus.params
#  decrease the migration rate in scenario 3...
scenario.list[[3]]@migration[[1]] <- scenario.list[[3]]@migration[[1]] * 0.1

# load scenarios
test.params@scenarios <- scenario.list

# set fastsimcoal check
test.params@sim.check.func <- fsc.scenarioCheck

# ---- Set analysis function ----
test.params@rep.analysis.func <- skeleSim::analysis_funcs

# ---- Run replicates ----
test.params <- runSim(test.params)
