# write fastsimcoal input file
fsc.write <- function(num.pops, Ne, sample.size = NULL, sample.times = NULL,
                      growth.rate = NULL, mig.rates = NULL,
                      hist.ev = NULL, num.chrom = 1,
                      locus.params = NULL, label = NULL) {

  if(is.null(label)) label <- "fastsimcoal.skeleSim"
  file <- paste(label, ".par", sep = "")
  mig.rates <- if(is.list(mig.rates)) mig.rates else list(mig.rates)
  hist.ev <- if(is.list(hist.ev)) do.call(rbind, hist.ev) else rbind(hist.ev)
  if (dim(hist.ev)[2]!=7) {hist.ev <- t(hist.ev)}   #somehow histev is becoming transposed...
  locus.params <- if(is.list(locus.params)) do.call(rbind, locus.params) else rbind(locus.params)
  if(nrow(locus.params) == 1 & num.chrom > 1) {
    locus.params <- do.call(rbind, lapply(1:num.chrom, function(i) locus.params[1, ]))
  }

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

  write("//Number of independent loci [chromosome]", file, append = T)
  write(paste(num.chrom, "1"), file, append = T)
  for(i in 1:num.chrom) {
    write("//Per chromosome: Number of linkage blocks", file, append = T)
    write("1", file, append = T)
    write("//Per block: data type, num loci, rec. rate and mut rate + optional parameters", file, append = T)
    write(paste(locus.params[i, ], collapse = " "), file, append = T)
  }

  invisible(file)
}
