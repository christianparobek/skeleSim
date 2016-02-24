#' @title Run fastsimcoal
#' @description Run fastsimcoal
#'
#' @param params a \linkS4class{skeleSim.params} object.
#'
#' @return a modified \linkS4class{skeleSim.params} object with the results of
#'   a fastsimcoal run.
#'
#' @export
#'
fsc.run <- function(params) {
  label <- currentLabel(params)
  sc <- currentScenario(params)

  # Check that folder is empty
  if(file.exists(params@wd)) for(f in dir(label, full.names = T)) file.remove(f)

  if(is.null(ploidy)) {
    pl <- attr(locus.params, "ploidy")
    ploidy <- if(is.null(pl)) 1 else pl
  }

  # Write fastsimcoal input file
  file <- fsc.write(
    num.pops = num.pops,
    Ne = pop.size,
    sample.size = sample.size,
    sample.times = sample.times,
    growth.rate = growth.rate,
    mig.rates = mig.rates,
    hist.ev = hist.ev,
    num.chrom = num.chrom,
    locus.params = locus.params,
    ploidy = ploidy,
    label = label
  )

  # Run fastsimcoal
  cmd.line <- paste(
    sc@simulator.params@fastsimcoal.exec, "-i", file, "-n 1",
    ifelse(params@quiet, "-q", ""), "-S -c0"
  )
  err <- if(.Platform$OS.type == "unix") {
    system(cmd.line, intern = F)
  } else {
    shell(cmd.line, intern = F)
  }

  if(err == 0) {
    if(!params@quiet) cat("fastsimcoal exited normally\n")
  } else {
    stop("fastsimcoal exited with error ", err, "\n")
  }

  chrom.pos <- if(is.null(num.chrom)) {
    tapply(locus.params$num.markers, locus.params$chromosome, sum)
  } else {
    rep(sum(locus.params$num.markers), num.chrom)
  }
  chrom.pos <- cumsum(chrom.pos)
  chrom.pos <- cbind(start = c(1, chrom.pos[-length(chrom.pos)] + 1), end = chrom.pos)

  arp.file <- file.path(label, paste(label, "_1_1.arp", sep = ""))
  params@rep.sample <- fsc.read(arp.file, chrom.pos, ploidy)
  params
}
