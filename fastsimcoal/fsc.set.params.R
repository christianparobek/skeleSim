fsc.histEvMat <- function(num.events = 1) {
  # -- historical events --
  # 1) Number of generations, t, before present at which the historical event
  #    happened
  # 2) Source deme (the first listed deme has index 0)
  # 3) Sink deme
  # 4) Expected proportion of migrants to move from source to sink.
  # 5) New size for the sink deme, relative to its size at generation t
  # 6) New growth rate for the sink deme
  # 7) New migration matrix to be used further back in time
  hist.ev <- c(
    num.gen = 10, source.deme = 1, sink.deme = 0, prop.migrants = 1,
    new.sink.size = 1, new.sink.growth = 0, new.mig.mat = 0
  )
  do.call(rbind, lapply(1:num.events, function(x) hist.ev))
}

fsc.locusParamsMat <- function(params) {
  locus.type <- switch(
    params@locus.type,
    microsat = "MICROSAT", snp = "SNP", sequence = "DNA"
  )
  if(locus.type == "DNA") params@num.loci <- 1
  num.loci <- params@num.loci

  # locus.length
  #  DNA: sequence length
  #  SNP or MICROSAT: number of loci
  locus.length <- params@sequence.length
  if(locus.type != "DNA") locus.length <- 1

  # mut.rate
  #   DNA: mutation rate per bp
  #   MICROSAT: mutation rate per locus
  #   SNP: minimum frequency for the derived allele
  mut.rate <- params@mut.rate

  # locus.param.5
  #   DNA: transition rate (1 / 3 = no bias)
  #   MICROSAT: geometric parameter for GSM (0 = SMM)
  locus.param.5 <- 1/3
  if(locus.type == "MICROSAT") locus.param.5 <- 0
  if(locus.type == "SNP") locus.param.5 <- NULL

  # locus.param.6: Number of different alleles for MICROSAT
  #   (0 = no range constraint)
  locus.param.6 <- 0
  if(locus.type %in% c("SNP", "DNA")) locus.param.6 <- NULL

  locus.type <- rep(locus.type, num.loci)

  # -- locus params --
  cbind(
    locus.type = locus.type,
    locus.length = locus.length,
    recomb.rate = 0,
    mut.rate = mut.rate,
    locus.param.5 = locus.param.5,
    locus.param.6 = locus.param.6
  )
}