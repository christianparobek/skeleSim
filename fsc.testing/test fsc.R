rm(list = ls())
source("fsc.classes.r")
source("fsc.histEvMat.r")
source("fsc.locusParams.r")
source("fsc.read.r")
source("fsc.write.r")

fsc.run <- function(num.pops, pop.size, sample.size, sample.times, growth.rate,
                    mig.rates, hist.ev, num.chrom, locus.params,
                    quiet = TRUE, fastsimcoal.exec = "fsc252", ploidy = NULL, label = NULL) {

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

  err <- system(paste(
    fastsimcoal.exec, "-i", file, "-n 1",
    ifelse(quiet, "-q", ""),
    "-S -c0"
  ), intern = F)

  chrom.pos <- if(is.null(num.chrom)) {
    tapply(locus.params$num.markers, locus.params$chromosome, sum)
  } else {
    rep(sum(locus.params$num.markers), num.chrom)
  }
  chrom.pos <- cumsum(chrom.pos)
  chrom.pos <- cbind(start = c(1, chrom.pos[-length(chrom.pos)] + 1), end = chrom.pos)

  label <- gsub(".par", "", file)
  arp.file <- file.path(label, paste(label, "_1_1.arp", sep = ""))
  fsc.read(arp.file, chrom.pos, ploidy)
}


num.pops <- 3
pop.size <- c(50, 100, 500)
sample.size <- c(25, 50, 25)
sample.times <- NULL
growth.rate <- NULL
mig.rates <- list(matrix(
  c(0, 0.01, 0.05, 0.025, 0, 0.025, 0.05, 0.01, 0),
  nrow = num.pops
))
hist.ev <- fsc.histEvMat(0)

dna.df <- do.call(rbind, lapply(1:8, function(x) {
  fsc.locus.dna(sample(5:30, 1), runif(1, 1e-8, 1e-2))
}))

msat.df <- do.call(rbind, lapply(1:4, function(x) {
  fsc.locus.msat(sample(5:20, 1), sample(5:10, 1), runif(1, 1e-10, 1e-3))
}))

snp.df <- do.call(rbind, lapply(1:6, function(x) {
  fsc.locus.snp(sample(1:4, 1), runif(1, 0, 0.01))
}))


num.chrom <- 3
locus.params <- snp.df

locus.params$chromosome <- sample(1:num.chrom, nrow(locus.params), T)
num.chrom <- NULL

result <- fsc.run(num.pops, pop.size, sample.size, sample.times, growth.rate,
                  mig.rates, hist.ev, num.chrom, locus.params)

print(result)
