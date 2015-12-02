#' @title Check all simulation parameters
#' @description Check all simulation parameters
#'
#' @param params a \linkS4class{skeleSim.params} object.
#'
#' @export
#'
overall.check <- function(params) {
  params@other.checks <- non.scenario.check(params)
  print(params@other.checks)
  #here we call the scenario checks (simulator specific and general)
  prv_chk<-params@sim.scen.checks  #store what is was in check slot
  #then calculate new checks
  ths_chk <- rbind(params@sim.check.func(params), gen.scenario.check(params))
  print(prv_chk);  print(ths_chk)
  #if what was there is null, replace with new checks
  if (is.null(prv_chk)) params@sim.scen.checks <- ths_chk
  #else, check which lines are there and replace info
  else {
    for (i in rownames(ths_chk)) {
      if (i %in% rownames(prv_chk)) prv_chk[i,] <- ths_chk[i,] #if it is there, replace it
      else {
        rbind(prv_chk,ths_chk[i,])  #if not, bind it
        rownames(prv_chk)[nrow(prv_chk)]<-i

      }
    }
    params@sim.scen.checks <- prv_chk
  }

  ##############TO DO write to a file error log#################

  #output result based on both sets of checks
  return(params)

}
