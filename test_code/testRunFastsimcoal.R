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

# create scenarios
scenario.list <- list(
  fsc.loadScenario(
    num.pops = 3,
    pop.size = c(50, 100, 500),
    sample.size = c(25, 50, 10),
    migration = list(matrix(
      c(0, 0.01, 0.05, 0.025, 0, 0.025, 0.05, 0.01, 0), nrow = 3
    )),
    locus.type = "dna",
    sequence.length = c(400, 100),
    mut.rate = c(1e-7, 1e-3)
  ),

  fsc.loadScenario(
    num.pops = 3,
    pop.size = c(50, 100, 500),
    sample.size = c(25, 50, 10),
    migration = list(matrix(
      c(0, 0.01, 0.05, 0.025, 0, 0.025, 0.05, 0.01, 0), nrow = 3
    )),
    locus.type = "msat",
    num.loci = c(10, 20, 10),
    mut.rate = c(1e-2, 1e-5, 1e-2),
    range.constraint = c(0, 10, 10)
  )
)

# load scenarios
test.params@scenarios <- scenario.list

# set fastsimcoal check
test.params@sim.check.func <- fsc.scenarioCheck

# ---- Set analysis function ----
test.params@rep.analysis.func <- skeleSim::analysis_funcs

# ---- Run replicates ----
test.params <- runSim(test.params)
