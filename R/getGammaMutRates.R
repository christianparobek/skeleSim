#' @name getGammaMutRates
#' @title Get a Vector of Randomly Chosen Mutation Rates
#' @description Get a vector of mutation rates from a gamma distribution by
#'   specifying the mean and std. dev. of the distribution
#'
#' @param n Number of loci to select rates for
#' @param gmean Mean of the gamma distribution
#' @param gstd Std Dev of the gamma distribution
#'
#' @return a vector of mutation rates (type numeric)
#'
#' @examples
#' rates = getGammaMutRates(1000,gmean=0.0001, gstd=0.0001)
#' hist(rates)
#' rates = getGammaMutRates(1000,gmean=0.0001, gstd=0.00001)
#' hist(rates)
#' @importFrom stats rgamma
#' @export
#'
getGammaMutRates <- function(n, gmean = 0.0001, gstd = 0.00001) {
    mu.mean <- gmean
    mu.sd <- gstd
    scale <- (mu.sd ^ 2) / mu.mean
    shape <- (mu.mean / mu.sd) ^ 2
    rgamma(n, scale = scale, shape = shape)
}
