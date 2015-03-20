stopifnot(require("reshape2"))
stopifnot(require("ggplot2"))
stopifnot(require("gridExtra"))

plot_all_stats<-function(main_list){
  
  #in the future this will be the real results dataframe
  results_df<-main_list$results_from_analyses
  num_scen<-length(main_list$scenarios_list[,1])
  stats_names<-names(main_list$results_from_analyses)
    
#   but for now we will make one up, assuming TWO scenarios
#   num_scen<-2
#   results_df<-data.frame(matrix(NA,nrow=main_list$common_params$num_reps*num_scen,ncol=3))
#   here I will put the names of the functions/summ statistics that get called
#   stats_names<-c("alleles","het","fis","gst")
#   
#   these results are totally synthetic
#   results_df[,1]<-c(rep(1,main_list$common_params$num_reps),rep(2,main_list$common_params$num_reps))
#   results_df[,2]<-abs(floor(rnorm(main_list$common_params$num_reps*num_scen,1,10)))
#   results_df[,3]<-runif(main_list$common_params$num_reps*num_scen)
#   results_df[,4]<-rnorm(main_list$common_params$num_reps*num_scen)
  colnames(results_df)<-c("scenario",stats_names)
  
  #plotting option 1
  results_melted<-melt(results_df,measure.vars=c(stats_names))
  ggplot(results_melted, aes(value)) + 
      geom_density() + facet_wrap(scenario ~ variable, 
      ncol = length(stats_names), scales = "free")
  
  #plotting option 2
  plots <- tapply(1:nrow(results_melted), list(results_melted$scenario, results_melted$variable), function(i) {
    df <- results_melted[i, ]
    sc <- unique(df$scenario)
    vr <- unique(df$variable)
    ggplot(df, aes(value)) + geom_density() + geom_vline(xintercept=table_means[sc, vr])
  })
  do.call(grid.arrange, plots)

}