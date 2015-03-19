#' @title Run FASTSIMCOAL
#' @description Run FASTSIMCOAL to generate a list of simulated gtypes.
#'
#' @param num.pops number of populations.
#' @param Ne effective population size.
#' @param sample.size number of samples to take.
#' @param sample.time time to draw samples.
#' @param growth.rate growth rate of populations.
#' @param mig.mat migration matrix.
#' @param hist.ev historical events.
#' @param num.chrom number of chromosomes.
#' @param data.type type of data.
#' @param locus.params locus parameters.
#' @param label character string to label files with.
#' @param num.sims number of simulations to run.
#' @param inf.site.model logical. Infinite site model?
#' @param quiet logical. Run quietly?
#'
#' @return Invisibly returns he name of the fastsimcoal .par file.
#'
#' @note Assumes that the program \code{fastsimcoal} is properly installed and
#'   available on the command line. On PC's, this requires having it in a folder in
#'   the PATH environmental variable. On Macs, the executable should be installed
#'   in a folder like \code{/usr/local/bin}.
#'
#' @references Excoffier, L. and Foll, M (2011) fastsimcoal: a continuous-time coalescent
#'   simulator of genomic diversity under arbitrarily complex evolutionary scenarios
#'   Bioinformatics 27: 1332-1334.  \url{http://cmpg.unibe.ch/software/fastsimcoal2/}
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export

fastsimcoal.skeleSim.run <- function(num.pops, Ne, sample.size = NULL, sample.time = NULL,
                            growth.rate = NULL, mig.rates = NULL,
                            hist.ev = NULL, num.chrom = 1,
                            locus.params = NULL, label = "fastsimcoal.skeleSim",
                            num.sims = 1, inf.site.model = TRUE, quiet = TRUE) {

  # Write input file
  hist.ev <- if(is.list(hist.ev)) do.call(rbind, hist.ev) else rbind(hist.ev)
  locus.params <- if(is.list(locus.params)) do.call(rbind, locus.params) else rbind(locus.params)
  if(nrow(locus.params) == 1 & num.chrom > 1) locus.params <- do.call(rbind, lapply(1:num.chrom, function(i) locus.params[1, ]))

  file <- paste(label, ".par", sep = "")
  write(paste("//  <<", label, ">>  (input from 'fastsimcoal.skeleSim.run')"), file)
  write(paste(num.pops, "populations to sample"), file, append = T)

  write("//Population effective sizes", file, append = T)
  for(i in 1:length(Ne)) write(Ne[i], file, append = T)

  write("//Sample sizes", file, append = T)
  if(is.null(sample.size)) sample.size <- rep(Ne, length(num.pops))
  if(is.null(sample.time)) sample.time <- rep(0, length(num.pops))
  sample.size <- paste(sample.size, sample.time)
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

  write("//Number of independent loci [chromosome]", file, append = T)
  write(paste(num.chrom, "0"), file, append = T)
  for(i in 1:num.chrom) {
    write("//Per chromosome: Number of linkage blocks", file, append = T)
    write("1", file, append = T)
    write("//Per block: data type, num loci, rec. rate and mut rate + optional parameters", file, append = T)
    write(paste(locus.params[i, ], collapse = " "), file, append = T)
  }

  # Check/setup folder structure
  if(file.exists(label)) for(f in dir(label, full.names = T)) file.remove(f)

  # Run fastsimcoal
  cmd <- paste("fastsimcoal -i", file, "-n", num.sims,
    ifelse(inf.site.model, "-I", ""), ifelse(quiet, "-q", ""), sep = " "
  )
  err <- system(cmd, intern = F)

  if(err == 0) {
    if(!quiet) cat("fastsimcoal exited normally\n")
  } else {
    stop("fastsimcoal exited with error ", err, "\n")
  }

  invisible(file)
}
