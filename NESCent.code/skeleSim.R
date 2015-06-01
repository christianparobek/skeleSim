rm(list = ls())
gc()
source("fastsimcoal.skeleSim.R")
source("genind.metadata.getter.R")
source("new.mainparam.list.R")
source("plot_all_stats.R")
source("summary_table_stats.R")

sim.choice <- function() {
  # question text
  questions <- c(
    snps = "Do you have SNP data?",
    non.diploid = "Is your data other than diploid?",
    marker.num = "Do you want to simulate many markers?",
    pop.size = "Do you have large population sizes?",
    complex.hist = "Do you have a complex history to simulate?",
    deep.time = "Are you looking at deep time frames?",
    demography = "Do you want to include demography?",
    management = "Does your question involve management decisions?",
    completion.time = "Do you need a short completion time?",
    computer = "Do you have large computer capacity?"
  )
  # default responses
  responses <- c(F, F, F, F, F, F, F, F, F, F)
  # response weights
  forward.wts <- c(0, 0, 0.3, 0.2, 0, 0.2, 1, 1, 0.2, 0.3)
  names(responses) <- names(forward.wts) <- names(questions)

  # loop through each question
  for(this.quest in names(questions)) {
    # create prompt with default
    prompt.def <- if(responses[this.quest]) "Y/n" else "y/N"
    prompt <- paste(questions[this.quest], " (", prompt.def, ") ", sep = "")
    # keep asking question until appropriate response received
    ans <- "empty"
    while(!ans %in% c("y", "n", "")) {
      ans <- tolower(readline(prompt))
    }
    # if not default, change response
    if(ans != "") responses[this.quest] <- ans == "y"
  }
  cat("\n")

  # find reasons that forward-time models are excluded
  fwd.excl <- forward.wts == 0 & responses
  if(any(fwd.excl)) {
    reasons <- paste(names(questions)[fwd.excl], collapse = ", ")
    cat("Forward-time simulations are excluded because: ", reasons, "\n")
  }

  # find reasons that forward-time models are required
  fwd.req <- forward.wts == 1 & responses
  if(any(fwd.req)) {
    reasons <- paste(names(questions)[fwd.req], collapse = ", ")
    cat("Forward-time simulations are required because: ", reasons, "\n")
  }

  # get relative 'score' for each model
  fwd.score <- sum(forward.wts * responses) / length(responses)
  cat("Coalescent score: ", 1 - fwd.score, "\n")
  cat("Forward-time score: ", fwd.score, "\n")

  ans <- "empty"
  prompt <- "Choose a simulator: (c)oalescent or (f)orward-time: "
  while(!ans %in% c("c", "f")) {
    ans <- tolower(readline(prompt))
  }

  return(ans)
}

# create multiple prompts to generate a vector
vec.prompt <- function(prompt, n) {
  sapply(1:n, function(i) readline(paste(prompt, " #", i, ": ", sep = "")))
}

# set migration rates
create.mig.matrix <- function(params) {
  stopifnot(require(rmetasim))

  params$common_params$overall_mig_rate <- as.numeric(readline("Overall migration rate: "))
  se <- tolower(readline("Do you want to define a spatially-explicit model? (y/N): ")) == "y"
  R.int <- if(se) {
    h.dim <- as.integer(readline("  Enter number of rows for landscape matrix: "))
    params$common_params$h.dim <- c(h.dim, num.pops / h.dim)
    if(!identical(trunc(h.dim[2]), h.dim[2])) {
      err.txt <- paste("You can't form a landscape matrix for", num_pops,
                       "populations with", h.dim[1], "rows.")
      stop(err.txt)
    }
    params$common_params$mean_mig_dist <- readline("  Enter mean migration distance: ")
    params$common_params$mig.model <- "distance"
    landscape.mig.matrix(h = params$common_params$num_pops,
                         h.dim = params$common_params$h.dim,
                         mig.model = params$common_params$mig.model,
                         distance.fun = dexp,
                         rate = 1 / params$common_params$mean_mig_dist,
                         distance.factor = 1
    )$R.int
  } else {
    params$common_params$mig.model <- "island"
    landscape.mig.matrix(h = params$common_params$num_pops,
                         mig.model = params$common_params$mig.model)$R.int
  }
  params$common_params$mig_rates <- R.int * params$common_params$overall_mig_rate

  return(params)
}

set.common.params <- function(params, genind.obj = NULL) {
  cat("\n")
  params$proj_title <- readline("Enter a title for the simulation: ")
  cat("\n")
  cat("--- Population Information ---\n")
  if(params$user_has_data) {
    params$common_params$num_pops <- genind.metadata.getter(gen_ind_obj)$NumberOfPops
    params$common_params$sample_sizes <- genind.metadata.getter(gen_ind_obj)$SampsPerPop
    cat("Number of populations: ", params$common_params$num_pops, "\n")
    cat("Number of samples from each population: ", params$common_params$sample.sizes, "\n")
  } else {
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
  }
  params <- create.mig.matrix(params)

  cat("\n")
  cat("--- Locus Information ---\n")
  locus.type <- ""
  prompt <- "Locus type (m)icrosatellite, (d)na sequence, (s)np: "
  while(!locus.type %in% c("m", "d", "s")) locus.type <- tolower(readline(prompt))
  locus.type <- switch(locus.type, m = "microsat", d = "sequence", s = "snp")
  params$common_params$locus_type <- locus.type
  if(params$user_has_data) {
    params$common_params$num_loci <- genind.metadata.getter(gen_ind_obj)$NumberOfLoci
    cat("Number of loci: ", params$common_params$num_loci, "\n")
  } else {
    num.loci <- as.integer(if(locus.type != "sequence") readline("Number of loci: ") else 1)
    params$common_params$num_loci <- num.loci
  }
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

  return(params)
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

run.sims <- function(params) {
  params <- sim.iterator(params)
  save(params, file = paste(params$proj_title, ".last.rep.rdata", sep = ""))
  plot_all_stats(params)
  params <- summary_table_stats(params)
  cat("/n")
  print(params$summary_results_table)
  save(params, file = paste(params$proj_title, ".results.rdata", sep = ""))
  return(params)
}

skeleSim <- function(quiet = FALSE, genind.obj = NULL) {
  # set main and common parameters
  params <- new.mainparam.list()
  params$quiet <- quiet
  params$sim_chosen <- sim.choice()
  params$user_has_data <- !is.null(genind.obj)
  params <- set.common.params(params, genind.obj)

  # set scenarios
  # toy scenarios_list << REMOVE THIS >>
  params$scenarios_list <- data.frame(Ne = 1)

  # set simulator-specific parameters
  params <- switch(params$sim_chosen,
                   c = set.fastsimcoal.params(params),
                   f = set.specparams.rmetasim.R(params)
  )

  # set analysis functions
  # toy analysis.func << REMOVE THIS >>
  params$analysis.func <- function(params){
    params$rep.analysis <- c("stat1"=abs(floor(rnorm(1,1,10))),
                             "stat2"=runif(1),
                             "stat3"=rnorm(1),
                             "stat4"=rbinom(1,1,.5))
    return(params)
  }

  # run simulator
  params <- run.sims(params)

  invisible(params)
}
