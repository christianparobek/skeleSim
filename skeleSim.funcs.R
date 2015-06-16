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
	file <- paste(label, ".params.rdata", sep = "")
	if(!dir.exists(label)) dir.create(label)
	save(params, file = file.path(label, file))
        c(scenario = sc, params@rep.result)
      }))
    }))
  }, finally = setwd(wd))

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
