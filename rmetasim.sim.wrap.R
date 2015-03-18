#' @title rmetasim.sim.wrap
#' A wrapper for rmetasim simulations to be used in the skeleSim framework.
#' @parlist The main input parameter list for the simulation.  This must contain the
#' element common_params and the element spec_params_rmetasim
#'
rmetasim.sim.wrap <- function(parlist)
    {
        if (! is.list(parlist)) {stop("parlist must be a list object")}
        if (! ("common_params" %in% names(parlist))) {stop ("need a common_params list")}
        if (! ("spec_params_rmetasim" %in% names(parlist))) {stop ("need a spec_params_rmetasim list")}

        
              
    }

