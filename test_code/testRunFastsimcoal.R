rm(list = ls())
library(skeleSim)
set.seed(1)

# ---- Load parameters ---

# create new skeleSim parameters object
test.params <- new("skeleSim.params")
test.params@title <- "testRun"
test.params@date <- Sys.time()
test.params@quiet <- FALSE
test.params@question <- "n"
test.params@simulator.type <- "c"
test.params@simulator <- "fsc"
test.params@num.reps <- 3
test.params@num.perm.reps <- 1000
test.params@num.cores <- 1
test.params@sim.func <- fsc.run
test.params@wd <- "testRun.wd"

# create and load scenarios
test.params@scenarios <- list(
  fsc.loadScenario(
    num.pops = 3,
    pop.size = c(20000, 5000, 10000),
    sample.size = c(20, 20, 6),
    sample.times = c(0, 0, 1500),
#     migration = list(matrix(
#       c(0, 0.01, 0.05, 0.025, 0, 0.025, 0.05, 0.01, 0), nrow = 3
#     )),
    hist.ev = fsc.histEvMat(
      num.gen = c(2000, 2980, 3000, 15000),
      source.deme = c(1, 1, 1, 0),
      sink.deme = c(2, 1, 0, 2),
      prop.migrants = c(0.05, 0, 1, 1),
      new.sink.size = c(1, 0.04, 1, 3),
      new.sink.growth = 0,
      new.mig.mat = 0
    ),
    locus.type = "dna",
    sequence.length = c(100, 200, 200, 300, 500),
    mut.rate = c(1e-4, 1e-5, 1e-6, 1e-7, 1e-3),
    chromosome = c(1, 1, 2, 2, 3)
  )
)

# set fastsimcoal check
test.params@sim.check.func <- fsc.scenarioCheck

# ---- Set analysis function ----
test.params@rep.analysis.func <- skeleSim::analysisFunc

# ---- Run replicates ----
result <- runSim(test.params)

save(result, file = paste(result$params@title, "result.rdata"))
