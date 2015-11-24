non.scenario.check <-
function(params) {
  results.check <- c(
    title.not.null = !is.null(params@title),
    #check that number of reps is greater than 0
    at.least.1.rep = params@num.reps > 0
  )
  print(results.check)
  return(results.check)
}
