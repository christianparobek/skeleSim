#' @title Create fastsimcoal historical event matrices
#' @description Create fastsimcoal historical event matrices
#'
#' @param num.events number of historical events.
#'
#' @return a blank fastsimcoal historical event matrices that can be
#'   filled in later
#'
#' @export
#'
fsc.histEvMat <- function(num.events = 0) {
  # -- historical events --
  # 1) Number of generations, t, before present at which the historical event
  #    happened
  # 2) Source deme (the first listed deme has index 0)
  # 3) Sink deme
  # 4) Expected proportion of migrants to move from source to sink.
  # 5) New size for the sink deme, relative to its size at generation t
  # 6) New growth rate for the sink deme
  # 7) New migration matrix to be used further back in time
  if(num.events == 0) return(NULL)
  hist.ev <- c(
    num.gen = 0, source.deme = 0, sink.deme = 0, prop.migrants = 1,
    new.sink.size = 1, new.sink.growth = 0, new.mig.mat = 0
  )
  do.call(rbind, lapply(1:num.events, function(x) hist.ev))
}
