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
  
  results_df<-params@analysis.results
  num_scen<-length(params@num.reps)
  stats_names<-names(params@analysis.results)
  colnames(results_df)<-c(stats_names)
  num_stats<-length(stats_names)-1
  
  #table of means and SD
  table_means<-data.frame(matrix(NA,nrow=num_scen,ncol=num_stats))
  table_sd<-data.frame(matrix(NA,nrow=num_scen,ncol=num_stats))
  colnames(table_means)<-stats_names[-1]; colnames(table_sd)<-stats_names[-1]
  for (i in 1:num_stats)
    table_means[,i]<-tapply(results_df[,1+i],results_df$scenario,mean)
  for (i in 1:num_stats)
    table_sd[,i]<-tapply(results_df[,1+i],results_df$scenario,sd)
  results_list<-list()
  results_list$means<-table_means
  results_list$sd<-table_sd
  params@summary.results<-results_list
  
  return(params)
}
