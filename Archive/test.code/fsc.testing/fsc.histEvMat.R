#' @title Create fastsimcoal historical event matrices
#' @description Create fastsimcoal historical event matrices
#'
#' @param num.gen Number of generations, t, before present at which the
#'   historical event happened.
#' @param source.deme Source deme (the first listed deme has index 0)
#' @param sink.deme Sink deme
#' @param prop.migrants Expected proportion of migrants to move from source to sink.
#' @param new.sink.size New size for the sink deme, relative to its size at
#'   generation t.
#' @param new.sink.growth New growth rate for the sink deme.
#' @param new.mig.mat New migration matrix to be used further back in time.
#' @param hist.ev a matrix describing historical events, with one row per event.
#' @param pop.size numerical vector giving size of each population.
#' @param growth.rate numerical vector giving growth rate of each population.
#' @param num.mig.mats number of migration matrices.
#'
#' @return a blank fastsimcoal historical event matrices that can be
#'   filled in later
#'
#' @export
#'
fsc.histEvMat <- function(num.gen = 0, source.deme = 0, sink.deme = 0,
                          prop.migrants = 1, new.sink.size = 1,
                          new.sink.growth = 0, new.mig.mat = 0) {
  cbind(
    num.gen = num.gen, source.deme = source.deme, sink.deme = sink.deme,
    prop.migrants = prop.migrants, new.sink.size = new.sink.size,
    new.sink.growth = new.sink.growth, new.mig.mat = new.mig.mat
  )
}


#' @rdname fsc.histEvMat
#'
fsc.histEvConverges <- function(hist.ev, pop.size, growth.rate, num.mig.mats = NULL) {
  if(is.null(hist.ev)) return(TRUE)
  if(length(growth.rate) == 1) growth.rate <- rep(growth.rate, length(pop.size))
  if(length(pop.size) != length(growth.rate)) {
    cat("'pop.size' and 'growth.rate' vectors must be same size.\n")
    return(FALSE)
  }
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


#' @rdname fsc.histEvMat
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