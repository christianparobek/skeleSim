#' @title write fsc files
#' @export
fscWrite <- function(pop.info, locus.params, mig.rates = NULL, hist.ev = NULL, label = NULL) {
  opt <- options(scipen = 999)

  ploidy <- attr(locus.params, "ploidy")
  pop.info[, c("pop.size", "sample.size")] <- pop.info[, c("pop.size", "sample.size")] * ploidy
print("one")
  if(is.null(label)) label <- "fastsimcoal.output"
  file <- paste(label, ".par", sep = "")
  mig.rates <- if(!is.null(mig.rates)) if(is.list(mig.rates)) mig.rates else list(mig.rates)
  hist.ev <- if(is.list(hist.ev)) do.call(rbind, hist.ev) else rbind(hist.ev)
  print("two")
  print(file)
  # Write input file
  print(getwd())

  write(paste("//  <<", label, ">>  (input from 'fastsimcoal.skeleSim.run')"), file=file)
  write(paste(nrow(pop.info), "populations to sample"), file=file, append = T)
  
  write("//Population effective sizes", file=file, append = T)
  for(i in 1:nrow(pop.info)) write(pop.info[i, "pop.size"], file=file, append = T)
    print("three")
  
  write("//Sample sizes", file=file, append = T)
  for(i in 1:nrow(pop.info)) write(pop.info[i, c("sample.size", "sample.times")], file=file, append = T)

  print(pop.info)
  
  write("//Growth rates", file=file, append = T)
  for(i in 1:nrow(pop.info)) write(pop.info[i, "growth.rate"], file=file, append = T)

  print("four")

  write("//Number of migration matrices", file=file, append = T)
  write(length(mig.rates), file=file, append = T)
  if(!is.null(mig.rates)) {
    for(i in 1:length(mig.rates)) {
      write("//migration matrix", file=file, append = T)
      for(r in 1:nrow(mig.rates[[i]])) write(mig.rates[[i]][r, ], file=file, append = T)
    }
  }

print("five")
  
  write("//Historical events: time, source, sink, migrants, new size, growth rate, migr. matrix", file=file, append = T)
  write(ifelse(is.null(hist.ev), 0, nrow(hist.ev)), file=file, append = T)
  if(!is.null(hist.ev)) {
    for(i in 1:nrow(hist.ev)) {
      write(paste(hist.ev[i, ], collapse = " "), file=file, append = T)
    }
  }

  print("six")
  
  num.chrom <- attr(locus.params, "num.chrom")
  if(!is.null(num.chrom)) locus.params$chromosome <- 1
  locus.params <- split(locus.params, locus.params$chromosome)
  num.independent <- if(is.null(num.chrom)) length(locus.params) else num.chrom
  chrom.struct <- if(num.independent == 1 | !is.null(num.chrom)) 0 else 1
  write("//Number of independent loci [chromosomes]", file=file, append = T)
  write(paste(num.independent, chrom.struct), file=file, append = T)
  for(block in locus.params) {
    block$chromosome <- NULL
    write("//Per chromosome: Number of linkage blocks", file=file, append = T)
    write(nrow(block), file=file, append = T)
    write("//Per block: data type, num loci, rec. rate and mut rate + optional parameters", file=file, append = T)
    for(i in 1:nrow(block)) {
      line <- paste(block[i, ], collapse = " ")
      write(gsub(" NA", "", line), file=file, append = T)
    }
  }

  print("seven")
  
  options(opt)
  invisible(file)
}
