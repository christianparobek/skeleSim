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

tic <- function(gcFirst = TRUE, type = c("elapsed", "user.self", "sys.self"), off = FALSE) {
  if(off) {
    assign(".type", NULL, envir = baseenv())
    assign(".tic", NULL, envir = baseenv())
    invisible(NULL)
  }
  type <- match.arg(type)
  assign(".type", type, envir = baseenv())
  if(gcFirst) gc(FALSE)
  tic <- proc.time()[type]
  assign(".tic", tic, envir = baseenv())
  invisible(tic)
}

toc <- function(show = FALSE) {
  type <- get(".type", envir = baseenv())
  if(is.null(type)) invisible(NULL)
  toc <- proc.time()[type]
  tic <- get(".tic", envir = baseenv())
  if(is.null(tic)) invisible(NULL)
  elapsed <- toc - tic
  if(show) print(elapsed)
  invisible(elapsed)
}

oneRep <- function(params) {
  # runs one replicate of simulator and loads sample into params@rep.sample
  params <- params@sim.func(params)
  # analyzes params@rep.sample and loads results into params@rep.result
  params <- params@rep.analysis.func(params)
  label <- currentLabel(params)
  file <- paste(label, ".params.rdata", sep = "")
  if(!dir.exists(label)) dir.create(label)
  save(params, file = file.path(label, file))
  c(scenario = params@current.scenario, params@rep.result)
}

stopRunning <- function(params, elapsed) {
  num.sc <- length(params@scenarios)
  num.reps <- params@num.reps
  reps.complete <- ((params@current.scenario - 1) * num.reps) + params@current.replicate
  total.reps <- num.sc * num.reps
  pct.complete <- round(100 * reps.complete / total.reps, 1)
  seconds.left <- (total.reps - reps.complete) * elapsed / num.reps
  eta <- Sys.time() + as.difftime(seconds.left, units = "secs")
  time.left <- eta - Sys.time()

  cat("\n--- skeleSim ---\n")
  cat("In ", round(elapsed, 1), " seconds, the simulator is ",
      pct.complete, "% complete: ",
      params@current.replicate, "/", num.reps, " replicates in ",
      params@current.scenario, "/", num.sc, " scenarios.\n", sep = "")
  cat("At this rate it is estimated to complete in:",
      round(time.left, 2), units(time.left), "\n")
  cont <- readline("Press 'n' to stop or any other key to continue: ")
  tolower(cont == "n")
}

runSim <- function(params) {
  # Check/setup folder structure
  if(file.exists(test.params@wd)) {
    unlink(test.params@wd, recursive = TRUE, force = TRUE)
  }
  dir.create(test.params@wd)
  wd <- setwd(test.params@wd)

  params@start.time <- Sys.time()
  params@analysis.results <- NULL
  tryCatch({
    num.sc <- length(params@scenarios)
    num.reps <- params@num.reps
    if(!is.null(params@timing)) tic()
      quit <- FALSE
      # must use nested for loops in order to be able to exit iterations
      for(sc in 1:num.sc) {
        if(quit) break # leave scenario for loop if user has decided to quit
        params@current.scenario <- sc
        for(r in 1:num.reps) {
          # run one replicate of simulator
          params@current.replicate <- r
          rep.result <- oneRep(params)
          # add analysis results to master matrix
          params@analysis.results <- rbind(params@analysis.results, rep.result)
          # check timing
          if(!is.null(params@timing)) {
            elapsed <- toc()
            if(!is.null(elapsed) & elapsed > params@timing) {
              params@timing <- NULL
              tic(off = TRUE)
              if(stopRunning(params, elapsed)) {
                quit <- TRUE
                break
              }
            }
          }
        }
      }

#       analysis.list <- lapply(1:length(params@scenarios), function(sc) {
#         params@current.scenario <- sc
#         do.call(rbind, lapply(1:params@num.reps, function(r) {
#           params@current.replicate <- r
#           oneRep(params)
#         }))
#       })
#       params@analysis.results <- do.call(rbind, analysis.list)

  }, finally = setwd(wd))
  rownames(params@analysis.results) <- NULL
  params@end.time <- Sys.time()

  cat("\n--- skeleSim ---\n")
  cat("Start:", format(params@start.time), "\n")
  cat("End:", format(params@end.time), "\n")
  elapsed <- round(params@end.time - params@start.time, 2)
  cat("Elapsed:", elapsed, units(elapsed), "\n")
  params
}

summ.stats.table<-function(params){

  results.datafr<-params@analysis.results
  num.sc <- length(params@scenarios)
  names.stats<-colnames(params@analysis.results)
  colnames(results.datafr)<-c(names.stats)
  num.stats<-length(names.stats)-1

  #table of means and SD
  table.means<-data.frame(matrix(NA,nrow=num.sc,ncol=num.stats))
  table.sd<-data.frame(matrix(NA,nrow=num.sc,ncol=num.stats))
  colnames(table.means)<-names.stats[-1]; colnames(table.sd)<-names.stats[-1]
  for (i in 1:num.stats)
    table.means[,i]<-tapply(results.datafr[,1+i],results.datafr[,1],mean)
  for (i in 1:num.stats)
    table.sd[,i]<-tapply(results.datafr[,1+i],results.datafr[,1],sd)
  results.list<-list()
  results.list$means<-table.means
  results.list$sd<-table.sd
  params@summary.results<-results.list

  return(params)
}



plot.all.stats<-function(params){
  stopifnot(require("reshape2"))
  stopifnot(require("ggplot2"))
  stopifnot(require("gridExtra"))

  results.datafr<-as.data.frame(params@analysis.results) 
  num.sc <- length(params@scenarios)
  names.stats<-colnames(params@analysis.results)
  colnames(results.datafr)<-c(names.stats)
  num.stats<-length(names.stats)-1
  
  #plotting option 1
  results.melted<-melt(results.datafr,id="scenario",measure.vars=c(names.stats[-1]))
  ggplot(results.melted, aes(value)) + 
      geom_density() + facet_wrap(scenario~variable, 
      ncol = num.stats, scales = "free")

  #histogram plot
  ggplot(results.melted, aes(value)) + geom_histogram() +
  facet_grid(variable~scenario, scales="free")
  
 }
