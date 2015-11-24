#' @title Check fastsimcoal historical events convergence
#'
#' @param hist.ev a matrix describing historical events, with one row per event.
#' @param pop.size numerical vector giving size of each population.
#' @param growth.rate numerical vector giving growth rate of each population.
#' @param num.mig.mats number of migration matrices.
#'
#' @return a logical stating whether convergence can be achieved given pattern of
#'   historical events.
#'
fsc.histEvConverges <- function(hist.ev, pop.size, growth.rate, num.mig.mats = NULL) {
  if(is.null(hist.ev)) return(TRUE)
  if(length(growth.rate) == 1) growth.rate <- rep(growth.rate, length(pop.size))
  if(length(pop.size) != length(growth.rate)) {
    cat("'pop.size' and 'growth.rate' vectors must be same size.\n")
    return(FALSE)
  }
  if(!fsc.histEvCheck(hist.ev, length(pop.size), num.mig.mats)) return(FALSE)
  hist.ev <- hist.ev[order(hist.ev[, 1], hist.ev[, 2], hist.ev[, 3]), , drop = FALSE]
  print(hist.ev)
  for(i in 1:nrow(hist.ev)) {
    gen <- hist.ev[i, 1]
    from <- hist.ev[i, 2] + 1
    to <- hist.ev[i, 3] + 1
    mig <- hist.ev[i, 4]
    new.size <- hist.ev[i, 5]
    pop.size[from] <- pop.size[from] * (1 - mig)
    t <- ifelse(i == 1, gen, gen - hist.ev[i - 1, 1])
    pop.size[to] <- pop.size[to] * exp(growth.rate[to] * t)
    growth.rate[to] <- hist.ev[i, 6]
    pop.size[to] <- pop.size[to] * new.size
    print(pop.size)
  }
  sum(pop.size > 0) == 1
}
