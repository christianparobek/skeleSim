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
test.params@simulator <- "f"
test.params@num.reps <- 10
test.params@timing <- 2
test.params@sim.func <- rms.run
test.params@wd <- "testRun.wd"

# create a base scenario parameters object. It will be copied and modified
#   later for different scenarios
base.scenario <- new("scenario.params")
base.scenario@num.pops <- 2
base.scenario@pop.size <- c(400,400)
base.scenario@sample.size <- 30
base.scenario@migration <- matrix(c(0,.1,.1,0),byrow=T,nrow=2) #Hard coded
base.scenario@locus.type <- "microsat"
base.scenario@num.loci <- 10
base.scenario@mut.rate <- rep(1e-4,10) #NOTE: this is a mutation rate for each locus

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
rms.params@init.pop.sizes = c(100,100,100,100)	#NOTE: this is per stage, per habitat, so four integers needed
rms.params@num.gen = 5
rms.params@num.alleles = rep(5,10) #NOTE this is a number of alleles for each locus
rms.params@allele.freqs = rep(.2,5)

# load rmetasim params
base.scenario@simulator.params <- rms.params

# create list of scenarios and modify
scenario.list <- lapply(1:3, function(i) base.scenario)
#  decrease the mutation rate in scenario 2...
scenario.list[[2]]@mut.rate <- 1e-5
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

# ---- Summarize analysis results ----
test.params <- summ.stats.table(test.params)
