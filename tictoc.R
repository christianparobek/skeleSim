tic <- function(gcFirst = TRUE, type = c("elapsed", "user.self", "sys.self")) {
  type <- match.arg(type)
  assign(".type", type, envir = baseenv())
  if(gcFirst) gc(FALSE)
  tic <- proc.time()[type]
  assign(".tic", tic, envir = baseenv())
  invisible(tic)
}

toc <- function(show = FALSE) {
  type <- get(".type", envir = baseenv())
  toc <- proc.time()[type]
  tic <- get(".tic", envir = baseenv())
  elapsed <- toc - tic
  if(show) print(elapsed)
  invisible(elapsed)
}

numIters <- function(seconds, params) {



runs <- 1
tic()
repeat {
  for (i in 1:50) mad(stats::runif(500))
  runs <- runs + 1
  if(toc() > 5) break
}
runs