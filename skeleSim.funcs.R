currentLabel <- function(params) {
  label <- paste(
    params@title,
    params@current.scenario,
    params@current.replicate,
    sep = "."
  )
  gsub("[[:punct:]]", ".", label)
}

currentScenario <- function(params) {
  params@scenarios[[params@current.scenario]]
}

runSim <- function(params) {
  # Check/setup folder structure
  if(file.exists(test.params@wd)) {
    unlink(test.params@wd, recursive = TRUE, force = TRUE)
  }
  dir.create(test.params@wd)
  wd <- setwd(test.params@wd)

  tryCatch({
    num.sc <- length(params@scenarios)
    params@analysis.results <- do.call(rbind, lapply(1:num.sc, function(sc) {
      params@current.scenario <- sc
      do.call(rbind, lapply(1:params@num.reps, function(r) {
        params@current.replicate <- r
        # runs one replicate of simulator and loads sample into params@rep.sample
        params <- params@sim.func(params)
        # analyzes params@rep.sample and loads results into params@rep.result
        params <- params@rep.analysis.func(params)
        label <- currentLabel(params)
        tryCatch({
		file <- paste(label, ".params.rdata", sep = "")
		dir.create(label)
	        save(params, file = file.path(label, file))
        	c(scenario = sc, params@rep.result)
		})
      }))
    }))
  }, finally = setwd(wd))

  params
}
