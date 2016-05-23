fsc.histEvMat <- function(num.events = 0) {
  # -- historical events --
  # 1) Number of generations, t, before present at which the historical event
  #    happened
  # 2) Source deme (the first listed deme has index 0)
  # 3) Sink deme
  # 4) Expected proportion of migrants to move from source to sink.
  # 5) New size for the sink deme, relative to its size at generation t
  # 6) New growth rate for the sink deme
  # 7) New migration matrix to be used further back in time
  if(num.events == 0) return(NULL)
  hist.ev <- c(
    num.gen = 0, source.deme = 0, sink.deme = 0, prop.migrants = 1,
    new.sink.size = 1, new.sink.growth = 0, new.mig.mat = 0
  )
  do.call(rbind, lapply(1:num.events, function(x) hist.ev))
}

fsc.locusParamsMat <- function(sc) {
  locus.type <- switch(
    sc@locus.type,
    microsat = "MICROSAT", snp = "SNP", sequence = "DNA"
  )
  num.loci <- sc@num.loci

  # locus.length
  #  DNA: sequence length
  #  SNP or MICROSAT: number of loci
  locus.length <- sc@sequence.length

  # mut.rate
  #   DNA: mutation rate per bp
  #   MICROSAT: mutation rate per locus
  #   SNP: minimum frequency for the derived allele
  mut.rate <- sc@mut.rate

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

fsc.histEvCheck <- function(hist.ev, num.pops, num.mig.mats = NULL) {
                                                                     
  if ((!is.numeric(as.matrix(hist.ev))) & (!(is.matrix(hist.ev)|is.data.frame(hist.ev)))) {
    cat("'hist.ev' must be a numerical matrix or numerical dataframe.\n")
    return(FALSE)
  }
  if(ncol(hist.ev) != 7) {
    cat("'hist.ev' must have 7 columns.\n")
    return(FALSE)
  }
  if(any(hist.ev[, 2] >= num.pops | hist.ev[, 3] >= num.pops)) {
    cat("values in columns 2 and 3 (source and sink demes) cannot be greater than the number of populations.\n")
    return(FALSE)
  }
  if(any(hist.ev[, 4] < 0 | hist.ev[, 5] < 0)) {
    cat("values in columns 4 (proportion of migrants) and 5 (new size for sink deme) must be greater than 0.\n")
    return(FALSE)
  }
  if(!is.null(num.mig.mats)) if(any(hist.ev[, 6] > num.mig.mats)) {
    cat("values in column 6 (new migration matrix) cannot be greater than the number of migration matrices.\n")
    return(FALSE)
  }
  TRUE
}

fsc.histEvConverges <- function(hist.ev, pop.size, growth.rate, num.mig.mats = NULL) {
  if(is.null(hist.ev)) return(TRUE)
  if(length(growth.rate) == 1) growth.rate <- rep(growth.rate, length(pop.size))
  if(length(pop.size) != length(growth.rate)) {
    cat("'pop.size' and 'growth.rate' vectors must be same size.\n")
    return(FALSE)
  }
  if(!fsc.histEvCheck(hist.ev, length(pop.size), num.mig.mats)) return(FALSE)
  hist.ev <- hist.ev[order(hist.ev[, 1], hist.ev[, 2], hist.ev[, 3]), , drop = FALSE]
  print(hist.ev)
  for(i in 1:nrow(hist.ev)) {
    gen <- hist.ev[i, 1]
    from <- hist.ev[i, 2] + 1
    to <- hist.ev[i, 3] + 1
    mig <- hist.ev[i, 4]
    new.size <- hist.ev[i, 5]
    pop.size[from] <- pop.size[from] * (1 - mig)
    t <- ifelse(i == 1, gen, gen - hist.ev[i - 1, 1])
    pop.size[to] <- pop.size[to] * exp(growth.rate[to] * t)
    growth.rate[to] <- hist.ev[i, 6]
    pop.size[to] <- pop.size[to] * new.size
    print(pop.size)
  }
  sum(pop.size > 0) == 1
}
