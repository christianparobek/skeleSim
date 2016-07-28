#' @title Create summary statistics table
#' @description Create summary statistics table
#'
#' @param params a \linkS4class{skeleSim.params} object.
#'
#' @importFrom stats sd
#'
summ.stats.table <- function(params) {
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
    table.sd[,i]<-tapply(results.datafr[,1+i],results.datafr[,1],stats::sd)
  results.list<-list()
  results.list$means<-table.means
  results.list$sd<-table.sd
  params@summary.results<-results.list

  return(params)
}
