summary_table_stats<-function(main_list){
  
  #in the future this will be the real results dataframe
  results_df<-main_list$results_from_analyses
  num_scen<-length(main_list$scenarios_list[,1])
  stats_names<-names(main_list$results_from_analyses)
  colnames(results_df)<-c("scenario",stats_names)
  
  #table of means and SD
  table_means<-data.frame(matrix(NA,nrow=num_scen,ncol=length(stats_names)))
  table_sd<-data.frame(matrix(NA,nrow=num_scen,ncol=length(stats_names)))
  colnames(table_means)<-stats_names; colnames(table_sd)<-stats_names
  for (i in 1:length(stats_names))
    table_means[,i]<-tapply(results_df[,1+i],results_df$scenario,mean)
  for (i in 1:length(stats_names))
    table_sd[,i]<-tapply(results_df[,1+i],results_df$scenario,sd)
  results_list<-list()
  results_list$means<-table_means
  results_list$sd<-table_sd
  
  return(results_list)
}