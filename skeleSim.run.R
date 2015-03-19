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

params <- list(
  label = "fastsimcoal.test",
  quiet = FALSE,
  user_has_data = FALSE
)
params <- set.commonparams(params)

main_list$common_params$num_pops<-num_pops
main_list$common_params$sample_sizes<-rep(20,num_pops)
main_list$common_params$num_loci<-10
main_list$common_params$pop_sizes<-rep(1000,num_pops)
main_list$common_params$overall_mig_rate <- 0.01
main_list$common_params$mut_rate<-0.0005
main_list$common_params$sequence_length <- 400
# number of simulation replicates to run
main_list$common_params$num_reps<-100

params <- set.fastsimcoal.params(params)

params$