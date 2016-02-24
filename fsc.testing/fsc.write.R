#' @title Write fastsimcoal input files
#' @description Write fastsimcoal input files
#'
#' @param num.pops number of populations.
#' @param Ne numeric vector: size of each population.
#' @param sample.size numeric vector: number of samples to draw from each population.
#' @param sample.times numeric vector: time step samples are drawn from for each population.
#' @param growth.rate numeric vector: growth rate of each population.
#' @param mig.rates list of matrices: migration rates between populations.
#' @param hist.ev matrix: one row per historical event.
#' @param num.chrom number of unique chromosomes.
#' @param locus.params data.frame: one row per locus parameter.
#' @param label character: label for file naming.
#'
fsc.write <- function(num.pops, Ne, sample.size = NULL, sample.times = NULL,
                      growth.rate = NULL, mig.rates = NULL,
                      hist.ev = NULL, num.chrom = NULL,
                      locus.params = NULL, ploidy = NULL, label = NULL) {

  if(is.null(ploidy)) {
    pl <- attr(locus.params, "ploidy")
    ploidy <- if(is.null(pl)) 1 else pl
  }
  Ne <- Ne * ploidy
  if(!is.null(sample.size)) sample.size <- sample.size * ploidy

  if(is.null(label)) label <- "fastsimcoal.skeleSim"
  file <- paste(label, ".par", sep = "")
  mig.rates <- if(is.list(mig.rates)) mig.rates else list(mig.rates)
  hist.ev <- if(is.list(hist.ev)) do.call(rbind, hist.ev) else rbind(hist.ev)

  # Write input file
  write(paste("//  <<", label, ">>  (input from 'fastsimcoal.skeleSim.run')"), file)
  write(paste(num.pops, "populations to sample"), file, append = T)

  write("//Population effective sizes", file, append = T)
  for(i in 1:length(Ne)) write(Ne[i], file, append = T)

  write("//Sample sizes", file, append = T)
  if(is.null(sample.size)) sample.size <- rep(Ne, length(num.pops))
  if(is.null(sample.times)) sample.times <- rep(0, length(num.pops))
  sample.size <- paste(sample.size, sample.times)
  for(i in 1:length(sample.size)) write(sample.size[i], file, append = T)

  write("//Growth rates", file, append = T)
  if(is.null(growth.rate)) growth.rate <- rep(0, num.pops)
  for(i in 1:length(growth.rate)) write(growth.rate[i], file, append = T)

  write("//Number of migration matrices", file, append = T)
  write(length(mig.rates), file, append = T)
  if(!is.null(mig.rates)) {
    for(i in 1:length(mig.rates)) {
      write("//migration matrix", file, append = T)
      for(r in 1:nrow(mig.rates[[i]])) write(mig.rates[[i]][r, ], file, append = T)
    }
  }

  write("//Historical events: time, source, sink, migrants, new size, growth rate, migr. matrix", file, append = T)
  write(ifelse(is.null(hist.ev), 0, nrow(hist.ev)), file, append = T)
  if(!is.null(hist.ev)) {
    for(i in 1:nrow(hist.ev)) {
      write(paste(hist.ev[i, ], collapse = " "), file, append = T)
    }
  }

  if(!is.null(num.chrom)) locus.params$chromosome <- 1
  locus.params <- split(locus.params, locus.params$chromosome)
  num.independent <- if(is.null(num.chrom)) length(locus.params) else num.chrom
  chrom.struct <- if(num.independent == 1 | !is.null(num.chrom)) 0 else 1
  write("//Number of independent loci [chromosomes]", file, append = T)
  write(paste(num.independent, chrom.struct), file, append = T)
  for(block in locus.params) {
    block$chromosome <- NULL
    write("//Per chromosome: Number of linkage blocks", file, append = T)
    write(nrow(block), file, append = T)
    write("//Per block: data type, num loci, rec. rate and mut rate + optional parameters", file, append = T)
    for(i in 1:nrow(block)) {
      line <- paste(block[i, ], collapse = " ")
      write(gsub(" NA", "", line), file, append = T)
    }
  }

  invisible(file)
}
