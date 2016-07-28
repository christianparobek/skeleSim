#' @title Create genpop object from landscape
#' @description Create genpop object from landscape
#'
#' @param l a \code{rmetasim} landscape object.
#'
#' @importFrom adegenet genind2genpop
#'
landscape.make.genpop <- function(l)
    {
        genind2genpop(landscape.make.genind(l))
    }
