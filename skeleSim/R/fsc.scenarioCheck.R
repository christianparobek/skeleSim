#' @title Check parameters for fastsimcoal
#' @description Check parameters for fastsimcoal
#'
#' @param params a \linkS4class{skeleSim.params} object.
#' @param hist.ev a matrix describing historical events, with one row per event.
#' @param pop.size numerical vector giving size of each population.
#' @param growth.rate numerical vector giving growth rate of each population.
#' @param num.mig.mats number of migration matrices.
#'
#' @export
#'
fsc.scenarioCheck <- function(params) {
  # check that sample times and growth rates are of length number of populations, and
  #   that historical events matrix converges
  results <- sapply(params@scenarios, function(sc) {
    fsc.histEvCheck(
      hist.ev = sc@simulator.params@hist.ev,
      pop.size = sc@pop.size, #simulator.params@pop.info[, "pop.size"],
      growth.rate = sc@simulator.params@growth.rate,
      num.mig.mats = length(sc@migration)
    )
  })
  results <- rbind(results)
  rownames(results) <- "hist.ev.good"
  colnames(results) <- paste("scenario", 1:ncol(results), sep = ".")
  return(results)
}

#' @rdname fsc.scenarioCheck
#' @export
#'
fsc.histEvConverges <- function(hist.ev, pop.size, growth.rate, num.mig.mats = NULL) {
  if(is.null(hist.ev)) return(TRUE)
  growth.rate <- rep(growth.rate, length.out = length(pop.size))
  hist.ev <- hist.ev[order(hist.ev[, 1], hist.ev[, 2], hist.ev[, 3]), , drop = FALSE]
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
  }
  sum(pop.size > 0) == 1
}


#' @rdname fsc.scenarioCheck
#' @export
#'
fsc.histEvCheck <- function(hist.ev, pop.size, growth.rate, num.mig.mats = NULL) {
  if(is.null(hist.ev)) return(TRUE)
  num.pops <- length(pop.size)
  if ((!is.numeric(as.matrix(hist.ev))) & (!(is.matrix(hist.ev)|is.data.frame(hist.ev)))) {
    cat("'hist.ev' must be a numerical matrix or numerical dataframe.\n")
    return(FALSE)
  }
  if(ncol(hist.ev) != 7) {
    cat("'hist.ev' must have 7 columns.\n")
    return(FALSE)
  }
  if(any(hist.ev[, 2] > (num.pops - 1) | hist.ev[, 3] > (num.pops - 1))) {
    cat("values in columns 2 and 3 (source and sink demes) cannot be greater than the number of populations.\n")
    return(FALSE)
  }
  if(any(hist.ev[, 4] < 0 | hist.ev[, 5] < 0)) {
    cat("values in columns 4 (proportion of migrants) and 5 (new size for sink deme) must be greater than 0.\n")
    return(FALSE)
  }
  if(!is.null(num.mig.mats)) if(any(hist.ev[, 6] > num.mig.mats)) {
    cat("values in column 6 (new migration matrix) cannot be greater than the number of migration matrices.\n")
    return(FALSE)
  }
  if(!fsc.histEvConverges(hist.ev, pop.size, growth.rate, num.mig.mats)) {
    cat("historical event matrix will not converge.\n")
    return(FALSE)
  }
  TRUE
}
