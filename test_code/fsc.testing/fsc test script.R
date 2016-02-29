rm(list = ls())
source("fsc code.r")

pop.info <- fsc.popInfo(
  pop.size = c(20000, 5000, 10000),
  sample.size = c(20, 20, 6),
  sample.times = c(0, 0, 1500)
)

locus.params <- fsc.locusParams(
  locus.type = "snp",
  mut.rate = 1e-5,
  num.loci = 1000
)

f <- fsc.run(
  pop.info,
  locus.params,
  mig.rates = matrix(
    c(0, 0.01, 0.05, 0.025, 0, 0.025, 0.05, 0.01, 0), nrow = 3
  ),
  hist.ev = fsc.histEvMat(
    num.gen = c(2000, 2980, 3000, 15000),
    source.deme = c(1, 1, 1, 0),
    sink.deme = c(2, 1, 0, 2),
    prop.migrants = c(0.05, 0, 1, 1),
    new.sink.size = c(1, 0.04, 1, 3),
    new.sink.growth = 0,
    new.mig.mat = 0
  )
)
