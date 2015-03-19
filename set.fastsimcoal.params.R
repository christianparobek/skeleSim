set.fastsimcoal.params <- function(params) {
  params$fastsimcoal.params <- list()
  params$fastsimcoal.params$sample.times <- NULL
  params$fastsimcoal.params$growth.rate <- NULL
  params$fastsimcoal.params$inf.site.model <- TRUE
  params$fastsimcoal.params$mig.rates <- list(params$common_params$mig_rates)

  # -- historical events --
  # 1) Number of generations, t, before prestent at which the historical even happened
  # 2) Source deme (the first listed deme has index 0)
  # 3) Sink deme
  # 4) Expected proportion of migrants to move from source to sink.
  # 5) New size for the sink deme, relative to its size at generation t
  # 6) New growth rate for the sink deme
  # 7) New migration matrix to be used further back in time
  num.gen <- 10
  source.deme <- 1
  sink.deme <- 0
  prop.migrants <- 1
  new.sink.size <- 1
  new.sink.growth <- 0
  new.mig.mat <- 0
  params$fastsimcoal.params$hist.ev <- cbind(
    num.gen, source.deme, sink.deme, prop.migrants,
    new.sink.size, new.sink.growth, new.mig.mat
  )

  # -- locus params --
  num.loci <- params$common_params$num_loci
  locus.type <- switch(as.character(params$common_params$locus_type),
    microsat = "MICROSAT", snp = "SNP", sequence = "DNA"
  )

  # locus.length
  #  DNA: sequence length
  #  SNP or MICROSAT: number of loci
  locus.length <- params$common_params$sequence_length
  if(locus.type != "DNA") locus.length <- 1

  # mut.rate
  #   DNA: mutation rate per bp
  #   MICROSAT: mutation rate per locus
  #   SNP: minimum frequency for the derived allele
  mut.rate <- params$common_params$mut_rate

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
  params$fastsimcoal.params$locus.params <- cbind(
    locus.type, locus.length, 0, mut.rate,
    locus.param.5, locus.param.6
  )

  return(params)
}
