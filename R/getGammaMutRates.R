getGammaMutRates <- function(n) {
  mut.dist.good <- FALSE
  scale <- shape <- NA
  while(!mut.dist.good) {
    mu.mean <- as.numeric(readline("  Gamma mean: "))
    mu.sd <- as.numeric(readline("  Gamma standard deviation: "))
    scale <- (mu.sd ^ 2) / mu.mean
    shape <- (mu.mean / mu.sd) ^ 2
    xlim <- qgamma(c(0.0001, 0.9999), scale = scale, shape = shape)
    if(length(unique(xlim)) == 1) xlim <- c(0, 1)
    curve(dgamma(x, scale = scale, shape = shape),
          xlab = "mutation rate", ylab = "density",
          xlim = xlim
    )
    ans <- tolower(readline("  Accept gamma parameters? (y/n)"))
    mut.dist.good <- ans == "y"
  }
  rgamma(n, scale = scale, shape = shape)
}