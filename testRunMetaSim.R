rm(list = ls())
source("skeleSim.classes.R")
source("rmetasim.skeleSim.R")
source("skeleSim.funcs.R")

# ---- Load parameters ---

# create new skeleSim parameters object
test.params <- new("skeleSim.params")
test.params@title <- "testRun"
test.params@date <- Sys.time()
test.params@quiet <- FALSE	#HARD CODE, for now
test.params@question <- "n"
test.params@simulator.type <- "f"	#input from user, choice
test.params@simulator <- "rms"
test.params@num.reps <- 10	#input from user, only integer allowed
test.params@timing <- 10		#input from user (or hard code), only integer allowed
test.params@sim.func <- rms.run
test.params@wd <- "testRun.wd"

# create a base scenario parameters object. It will be copied and modified
#   later for different scenarios
base.scenario <- new("scenario.params")		
base.scenario@num.pops <- 2		#input from user, only integer allowed
base.scenario@pop.size <- c(400,400)	#input from user, only vector of integers allowed
base.scenario@sample.size <- c(30,30)		#input from user, only integer allowed; could change to vector of integers
base.scenario@migration <- list(
          matrix(c(0,.1,.1,0),
          byrow=T,nrow=2
          )) 	#input from user, only list of matrices
base.scenario@locus.type <- "microsat"	#input from user, a choice of three types 
base.scenario@num.loci <- 10		#input from user, only integer allowed
base.scenario@mut.rate <- rep(1e-4,10) 		#input from user, only numeric between 0 and 1 allowed
						#NOTE: this is a mutation rate for each locus

# create rmetasim params object to load into base scenario
rms.params <- new("rmetasim.params")
rms.params@num.stgs <- 2		#HARD CODE, for now
rms.params@selfing = 0			#input from user, only numeric between 0 and 1 allowed
		#the next three are input from user, only matrix
rms.params@surv.matr = matrix(c(	0.2,0.0,
 			                        0.2,0.7),nrow=2,byrow=T)
rms.params@repr.matr = matrix(c(	0,10,
		                        	0,0),nrow=2,byrow=T)
rms.params@male.matr = matrix(c(	0,0,
	                        		0,1),nrow=2,byrow=T)
rms.params@init.pop.sizes = c(100,100,100,100)	#input from user, only vector of integers
					#NOTE: this is per stage, per habitat, so SxH integers needed
rms.params@num.gen = 1				#input from user, only vector of integers
rms.params@num.alleles = rep(5,10) 	#input from user, only vector of length number of loci, integers
					#NOTE this is a number of alleles for each locus
rms.params@allele.freqs = rep(.2,5)	#input from user, only vector of length number of alleles, numeric
					#NOTE currently assumes all loci have same number of alleles
					#and allele frequencies

# load rmetasim params
base.scenario@simulator.params <- rms.params

# create list of scenarios and modify
# user will need to be presented with list of modifiable parameters, which they choose one
# then they will choose a multiplier?
# can they choose more than one parameter?
scenario.list <- lapply(1:3, function(i) base.scenario)
#  decrease the mutation rate in scenario 2...
scenario.list[[2]]@mut.rate <- rep(1e-5,1)
#  decrease the migration rate in scenario 3...
scenario.list[[3]]@migration[[1]] <- scenario.list[[3]]@migration[[1]] * 0.1

# load scenarios
test.params@scenarios <- scenario.list


# ---- Set analysis function ----
# choices from user..
test.params@rep.analysis.func <- function(params) {
  result = rnorm(5)
  names(result) <- paste("result", 1:5, sep = ".")
  params@rep.result <- result
  params
}

# --- Set parameter check function for specific simulator---
test.params@sim.check.func <- rms.scenarioCheck


# ---- Run replicates ----
test.params <- runSim(test.params)


# ---- Summarize analysis results ----
#test.params <- summ.stats.table(test.params)

#plot.all.stats(test.params)

