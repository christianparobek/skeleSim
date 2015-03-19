rm(list = ls())
gc()
source("fastsimcoal.skeleSim.run.r")
source("genind.metadata.getter.r")
source("set.commonparams.r")
source("set.fastsimcoal.params.r")
source("sim.choice.r")
source("sim.wrap.fastsimcoal.r")
source("fastsimcoal.skeleSim.read.r")

vec.prompt <- function(prompt, n) {
  sapply(1:n, function(i) readline(paste(prompt, " #", i, ": ", sep = "")))
}

skeleSim.run <- function() {
  sim.type <- sim.choice()
  cat("\n")
  label <- readline("Enter a label for the simulation: ")

  params <- list(
    sim.type = sim.type,
    label = label,
    quiet = FALSE,
    user_has_data = FALSE
  )
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

  if(sim.type == "c") {
    params <- set.fastsimcoal.params(params)
    params$common_params$sim.func <- sim.wrap.fastsimcoal
  } else {
    params <- set.specparams.rmetasim.R(params)
    params$common_params$sim.func <- sim.wrap.rmetasim
  }

  save(params, file = paste(params$label, ".params.rdata", sep = ""))
  return(params)
}