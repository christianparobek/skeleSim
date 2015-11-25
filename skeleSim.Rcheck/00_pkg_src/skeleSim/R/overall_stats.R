#' @title Compute overall statistics
#' @description Compute overall statistics
#'
#' @param results_gtype a \code{\link{gtypes}} object.
#'
#' @importFrom strataG overallTest
#'
overall_stats <- function(results_gtype) {
  ovl <- overallTest(results_gtype, nrep = 5, quietly = TRUE)
  ovl.result <- ovl$result[complete.cases(ovl$result[,1]),]

  pnam <- c()
  for(i in 1:nrow(ovl.result))
    pnam <- c(pnam,rownames(ovl.result)[i],paste(rownames(ovl.result)[i],"pval", sep = ""))

  global.wide <- as.vector(t(ovl.result))
  names(global.wide) <- pnam
  global.wide
}
