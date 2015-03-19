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

skeleSim.run <- function() {
  sim.type <- sim.choice()
  label <- readline("Enter a label for the simulation: ")

  params <- list(
    sim.type = sim.type,
    label = label,
    quiet = FALSE,
    user_has_data = FALSE
  )
  params <- set.commonparams(params)

  cat("\n\n")
  cat("--- Population Information ---\n")
  num.pops <- readline("Number of populations: ")
  params$common_params$num_pops <- num.pops
  params$common_params$pop_sizes <- rep(readline("Size of each population: "), num.pops)
  params$common_params$sample_sizes <- rep(readline("Number of samples from each population: "), num.pops)
  params$common_params$overall_mig_rate <- readline("Overall migration rate: ")
  cat("\n\n")
  cat("--- Locus Information ---\n")
  locus.type <- ""
  prompt <- "Locus type (m)icrosatellite, (d)na sequence, (s)np: "
  while(!locus.type %in% c("m", "d", "s")) {
    locus.type <- tolower(readline(prompt))
  }
  locus.type <- switch(locus.type, m = "microsat", d = "sequence", s = "snp")
  main_list$common_params$locus_type <- locus.type
  params$common_params$num_loci <- readline("Number of loci: ")
  params$common_params$mut_rate <- readline("Mutation rate: ")
  if(locus.type == "sequence") {
    params$common_params$sequence_length <- readline("Sequence length: ")
  }
  cat("\n\n")
  cat("--- Simulation Information ---\n")
  params$common_params$num_reps <- readline("Number of replicates to run: ")

  if(sim.type == "c") {
    params <- set.fastsimcoal.params(params)
    params$common_params$sim.func <- sim.wrap.fastsimcoal
  } else {
    params <- set.specparams.rmetasim.R(params)
    params$common_params$sim.func <- sim.wrap.rmetasim
  }

  save(params, file = paste(param$label, ".params.rdata", sep = ""))
}