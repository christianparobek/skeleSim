#' @title Check fastsimcoal historical events
#' @description Test structure and consistency of fastsimcoal historical events object.
#'
#' @param hist.ev a matrix describing historical events, with one row per event.
#' @param num.pops number of populations.
#' @param num.mig.mats number of migration matrices.
#'
fsc.histEvCheck <- function(hist.ev, num.pops, num.mig.mats = NULL) {
  if ((!is.numeric(as.matrix(hist.ev))) & (!(is.matrix(hist.ev)|is.data.frame(hist.ev)))) {
    cat("'hist.ev' must be a numerical matrix or numerical dataframe.\n")
    return(FALSE)
  }
  if(ncol(hist.ev) != 7) {
    cat("'hist.ev' must have 7 columns.\n")
    return(FALSE)
  }
  if(any(hist.ev[, 2] >= num.pops | hist.ev[, 3] >= num.pops)) {
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
  TRUE
}
