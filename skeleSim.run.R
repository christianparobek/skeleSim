rm(list = ls())
gc()
source("fastsimcoal.skeleSim.R")
source("genind.metadata.getter.R")
source("new.mainparam.list.R")
source("set.commonparams.R")
source("sim.choice.R")
source("plot_all_stats.R")
source("summary_table_stats.R")

vec.prompt <- function(prompt, n) {
  sapply(1:n, function(i) readline(paste(prompt, " #", i, ": ", sep = "")))
}

sim.iterator <- function(megalist){
  ## Number of Scenarios
  num_scenarios <- nrow(megalist$scenarios_list)
  ## Number of Reps
  num_reps <- megalist$common_params$num_reps
  ## Define a "results_from_analysis" list
  megalist$results_from_analyses <-
    as.data.frame(do.call(rbind, lapply(1:num_scenarios, function(scenario) {
      megalist$current_scenario <- scenario
      do.call(rbind, lapply(1:num_reps, function(rep) {
        megalist$current_replicate <- rep
        megalist <- megalist$sim.func(megalist)
        megalist <- megalist$analysis.func(megalist)
        c(scenario = scenario, megalist$rep.analysis)
      }))
    })))
  return(megalist)
}


skeleSim.run <- function(quiet = FALSE) {
  params <- new.mainparam.list()
  params$quiet <- quiet
  params$sim_chosen <- sim.choice()

  cat("\n")
  params$proj_title <- readline("Enter a title for the simulation: ")

  params <- set.commonparams(params)

  cat("\n")
  cat("--- Population Information ---\n")
  num.pops <- as.integer(readline("Number of populations: "))
  params$common_params$num_pops <- num.pops
  pop.size.good <- FALSE
  while(!pop.size.good) {
    pop.sizes <- readline("Size of each population (press Enter to enter sizes individually): ")
    params$common_params$pop_sizes <- if(pop.sizes == "") {
      as.integer(vec.prompt("  Enter size of population", num.pops))
    } else rep(as.integer(pop.sizes), num.pops)
    sample.sizes <- readline("Number of samples from each population (press Enter to enter samples individually): ")
    params$common_params$sample_sizes <- if(sample.sizes == "") {
      as.integer(vec.prompt("  Enter number of samples from population", num.pops))
    } else rep(as.integer(sample.sizes), num.pops)
    pop.size.good <- all(as.integer(pop.sizes) >= as.integer(sample.sizes))
    if(!pop.size.good) cat("<< Error: Some sample sizes are larger than population sizes. Please re-enter. >>\n")
  }
  params$common_params$overall_mig_rate <- as.numeric(readline("Overall migration rate: "))

  cat("\n")
  cat("--- Locus Information ---\n")
  locus.type <- ""
  prompt <- "Locus type (m)icrosatellite, (d)na sequence, (s)np: "
  while(!locus.type %in% c("m", "d", "s")) locus.type <- tolower(readline(prompt))
  locus.type <- switch(locus.type, m = "microsat", d = "sequence", s = "snp")
  params$common_params$locus_type <- locus.type
  num.loci <- as.integer(if(locus.type != "sequence") readline("Number of loci: ") else 1)
  params$common_params$num_loci <- num.loci
  mut.rate <- readline("Mutation rate (press Enter to enter parameters of a Gamma distribution): ")
  mut.rate <- if(mut.rate == "") {
    mut.dist.good <- FALSE
    scale <- shape <- NA
    while(!mut.dist.good) {
      mu.mean <- as.numeric(readline("  Gamma mean: "))
      mu.sd <- as.numeric(readline("  Gamma standard deviation: "))
      scale <- (mu.sd ^ 2) / mu.mean
      shape <- (mu.mean / mu.sd) ^ 2
      curve(dgamma(x, scale = scale, shape = shape),
            xlab = "mutation rate", ylab = "density",
            xlim = c(0, 0.1)
      )
      ans <- tolower(readline("  Accept gamma parameters? (y/n)"))
      mut.dist.good <- ans == "y"
    }
    rgamma(num.loci, scale = scale, shape = shape)
  } else mut.rate
  params$common_params$mut_rate <- as.numeric(mut.rate)
  if(locus.type == "sequence") {
    params$common_params$sequence_length <- as.integer(readline("Sequence length: "))
  }

  cat("\n")
  cat("--- Simulation Information ---\n")
  params$common_params$num_reps <- as.integer(readline("Number of replicates to run: "))

  if(params$sim_chosen == "c") {
    params <- set.fastsimcoal.params(params)
    params$sim.func <- sim.wrap.fastsimcoal
    #params$analysis.func <- test.func
  } else {
    params <- set.specparams.rmetasim.R(params)
    params$sim.func <- sim.wrap.rmetasim
    #params$analysis.func <- test.func
  }

  # test analysis.func
  params$analysis.func <- function(params){
    params$rep.analysis <- c("stat1"=abs(floor(rnorm(1,1,10))),
                  "stat2"=runif(1),
                  "stat3"=7.1,
                  "stat4"=0.3)
    return(params)
  }

  save(params, file = paste(params$proj_title, ".params.rdata", sep = ""))

  params$scenarios_list <- data.frame(Ne = 1)
  params <- sim.iterator(params)
  plot_all_stats(params)
  params <- summary_table_stats(params)
  cat("/n")
  print(params$summary_results_table)
  save(params, file = paste(params$proj_title, ".results.rdata", sep = ""))
  invisible(params)
}
