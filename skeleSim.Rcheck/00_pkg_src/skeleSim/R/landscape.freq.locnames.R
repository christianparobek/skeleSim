#' @title Landscape locus names
#' @description Landscape locus names
#'
#' @param l a \code{rmetasim} landscape object.
#'
#'@importFrom rmetasim landscape.ploidy landscape.locus landscape.democol
#'
landscape.freq.locnames <- function(l)
  {
    num.loc <- length(landscape.ploidy(l))
    namevec <- vector("character", num.loc)
    for (loc in 1:num.loc)
      {
        genos <- landscape.locus(loc,l)[,-1:-landscape.democol()]
        loc.names <- paste(loc, names(table(unlist(genos))), sep = ".")
        namevec[loc] <- paste("L", loc.names, sep = '')
      }
    namevec
  }
