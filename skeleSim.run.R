rm(list = ls())
gc()
source("fastsimcoal.skeleSim.run.r")
source("genind.metadata.getter.r")
source("set.commonparams.r")
source("set.fastsimcoal.params.r")
source("sim.choice.r")
source("sim.wrap.fastsimcoal.r")
source("fastsimcoal.skeleSim.read.r")
library(rmetasim)

sim.type <- sim.choice()
label <- readline("Enter a label for the simulation: ")

params <- list(
  sim.type = sim.type,
  label = label,
  quiet = FALSE,
  user_has_data = FALSE
)
params <- set.commonparams(params)

# number of populations
params$common_params$num_pops <- num_pops
# size of each population
params$common_params$pop_sizes <- rep(1000, num_pops)
# number of samples per population
params$common_params$sample_sizes <- rep(20, num_pops)
# migration rate
params$common_params$overall_mig_rate <- 0.01

# locus type
main_list$common_params$locus_type <- "microsat"
# number of loci
params$common_params$num_loci <- 10
# mutation rate
params$common_params$mut_rate <- 0.0005
# if DNA sequences, number of base pairs
params$common_params$sequence_length <- 400

# number of simulation replicates to run
params$common_params$num_reps <- 100

if(
params <- set.fastsimcoal.params(params)

params$