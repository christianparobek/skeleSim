#' @name skeleSim.internals
#' @title SkeleSim internal functions
#' @description SkeleSim internal functions
#'
#' @param params a \linkS4class{skeleSim.params} object.
#'
#' @return
#'   \tabular{ll}{
#'     \code{currentScenario} \tab the parameters for the current scenario.\cr
#'     \code{currentLabel} \tab a character label representing current scenario
#'       and replicate.\cr
#'   }
#'
currentScenario <- function(params) {
  params@scenarios[[params@current.scenario]]
}

#' @rdname skeleSim.internals
#'
currentLabel <- function(params) {
  label <- paste(
    params@title,
    params@current.scenario,
    params@current.replicate,
    sep = "."
  )
  gsub("[[:punct:]]", ".", label)
}

# #' @rdname skeleSim.internals
# #'
# results2gtypes <- function(params) {
#   #params@rep.sample is either a genind or a list of DNAbin objects
#   if(inherits(params@rep.sample, "genind")) {
#     return(genind2gtypes(params@rep.sample))
#   } else if(inherits(params@rep.sample, "gtypes")) {
#     return(params@rep.sample)
#   } else if(is.list(params@rep.sample)) {
#     # Convert the list of DNAbin objects to gtypes
#     genes <- params@rep.sample
#     g <- sequence2gtypes(genes$dna.seqs, strata = genes$strata)
#     return(labelHaplotypes(g)$gtype)
#   } else {
#     warning("cannot convert results in params to gtypes, returning NULL")
#   }
#   return(NULL)
# }