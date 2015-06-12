rm(list = ls())
source("skeleSim.classes.R")
source("rmetasim.skeleSim.R")
library("rmetasim")
source("skeleSim.funcs.R")

# ---- Load parameters ---

# create new skeleSim parameters object
test.params <- new("skeleSim.params")
test.params@title <- "testRun"
test.params@date <- Sys.time()
test.params@quiet <- FALSE
test.params@question <- "n"
test.params@simulator <- "c"
test.params@num.reps <- 10
test.params@sim.func <- rms.run
test.params@wd <- "testRun.wd"

# create a base scenario parameters object. It will be copied and modified
#   later for different scenarios
base.scenario <- new("scenario.params")
base.scenario@num.pops <- 3
base.scenario@pop.size <- c(50, 100, 500)
base.scenario@sample.size <- c(25, 50, 25)
base.scenario@migration <- matrix(
  c(0, 0.01, 0.05, 0.025, 0, 0.025, 0.05, 0.01, 0),
  nrow = base.scenario@num.pops
)
base.scenario@locus.type <- "microsat"
base.scenario@num.loci <- 10
base.scenario@mut.rate <- 1e-4

# create rmetasim params object to load into base scenario
rms.params <- new("rmetasim.params")
rms.params@num.stgs <- 2
rms.params@selfing = 0
rms.params@surv.matr = matrix(c(	0.2,0.0,
 			                        0.2,0.7),nrow=2,byrow=T)
rms.params@repr.matr = matrix(c(	0,10,
		                        	0,0),nrow=2,byrow=T)
rms.params@male.matr = matrix(c(	0,0,
	                        		0,1),nrow=2,byrow=T)
rms.params@init.pop.size = c(1000,1000)
rms.params@num.gen = 20
rms.params@num.alleles = 5
rms.params@allele.freqs = c(rep,.2,5)

# load rmetasim params
base.scenario@simulator.params <- rms.params

# create list of scenarios and modify
scenario.list <- lapply(1:3, function(i) base.scenario)
#  decrease the mutation rate in scenario 2...
scenario.list[[2]]@mut.rate <- 1e-5
scenario.list[[2]]@simulator.params@locus.params <- fsc.locusParamsMat(scenario.list[[2]])
#  decrease the migration rate in scenario 3...
scenario.list[[3]]@migration <- scenario.list[[3]]@migration * 0.1

# load scenarios
test.params@scenarios <- scenario.list


# ---- Set analysis function ----
test.params@rep.analysis.func <- function(params) {
  result = rnorm(5)
  names(result) <- paste("result", 1:5, sep = ".")
  params@rep.result <- result
  params
}


# ---- Run replicates ----
test.params <- runSim(test.params)
