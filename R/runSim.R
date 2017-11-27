#' @title Run simulation
#' @description Run simulation
#'
#' @param params a \linkS4class{skeleSim.params} object.
#' @param num.secs number of seconds to run timing checks
#'
#' @importFrom swfscMisc autoUnits
#' @importFrom utils write.table write.csv
#'
#' @export
#'
runSim <- function(params, num.secs = NULL) {
  # these are parameters that are going to be hosed anyway.
  params@analysis.results <- NULL
  params@sim.scen.checks <- NULL
  params@rep.result <- NULL

  # check parameters
  params <- overall.check(params)
  # print(params@other.checks)
  # print(params@sim.scen.checks)
  if(!all(params@other.checks, params@sim.scen.checks)) {
    filewr <- "error.log"
    write.csv(params@sim.scen.checks, file = filewr)
    #print(params@other.checks)
    write("\n", file = filewr, append = T)
    write.table(params@other.checks, file = filewr, append = T)
    stop("parameters do not pass checks; see error log for details")
  }
  # cat("\nparameter check complete\n\n")

  # Check/setup folder structure
  if(file.exists(params@wd)) unlink(params@wd, recursive = TRUE, force = TRUE)
  dir.create(params@wd)
  wd <- setwd(params@wd)

  results <- list(timing = list())
  params@analysis.results <- NULL
  params@analyses.requested <- analyses.check(params@analyses.requested)
  tryCatch({
    num.sc <- length(params@scenarios)
    num.reps <- params@num.sim.reps
    params@scenario.reps <- as.matrix(expand.grid(scenario = 1:num.sc, replicate = 1:num.reps))
    quit <- FALSE
    # loop through replicates for scenarios
    num.iter <- nrow(params@scenario.reps)
    results$timing$start.time <- Sys.time()
    for(i in 1:num.iter) {
        #if (debug)
      cat("Scenario replicate", i, "\n")
      if(quit) break # leave replicate for loop if user has decided to quit
      params@current.scenario <- params@scenario.reps[i, "scenario"]
      params@current.replicate <- params@scenario.reps[i, "replicate"]
      # run one replicate of simulator
      cat("  Starting sim function...\n")
      params <- params@sim.func(params)
      # analyzes params@rep.sample and loads results into params@rep.result
      params <- params@rep.analysis.func(params)
      cat("  Analysis done\n")
      # -->> REMOVE FOR RELEASE: SAVING params OBJECT FOR TESTING <<--
      # label <- currentLabel(params)
      # file <- paste(label, ".params.rdata", sep = "")
      # if(!dir.exists(label)) dir.create(label)
      # save(params, file = file.path(label, file))
      #-----
      # check timing
      results$timing$end.time <- Sys.time()
      if(!is.null(num.secs)) {
        elapsed <- results$timing$end.time - results$timing$start.time
        units(elapsed) <- "secs"
        if(elapsed > num.secs) break
      }
    }

    elapsed <- results$timing$end.time - results$timing$start.time
    results$timing$completion.time <- autoUnits(num.iter * elapsed / i)
    results$timing$pct.complete <- round(100 * i / num.iter, 1)
    results$params <- params
  }, finally = setwd(wd))

  results[[2]]@timing <- results[[1]]
  results[[2]]
}
