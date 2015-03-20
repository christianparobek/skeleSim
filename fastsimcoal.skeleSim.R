# read a .arp file written by fastsimcoal
# filename (and path) is in params$fastsimcoal.params$arp.file
fastsimcoal.skeleSim.read <- function(params) {
  stopifnot(require(ape))
  f <- readLines(params$fastsimcoal.params$arp.file)

  # get start and end points of data blocks
  start <- grep("SampleData=", f) + 1
  end <- which(f == "}") - 2
  pos <- cbind(start, end)
  # compile data for each population
  pop.data <- do.call(rbind, lapply(1:nrow(pos), function(i) {
    f.line <- f[pos[i, 1]:pos[i, 2]]
    f.line <- gsub("[[:space:]]+", "--", f.line)
    data.mat <- do.call(rbind, strsplit(f.line, "--"))[, -2]
    data.mat <- cbind(rep(paste("Sample", i), nrow(data.mat)), data.mat)
  }))

  # get data type
  data.type <- f[grep("DataType=", f)]
  data.type <- gsub("\tDataType=", "", data.type)
  is.seq <- switch(data.type, DNA = T, MICROSAT = F, STANDARD = F)
  result <- if(is.seq) {
    # replace sequence with all A's if there are no variable sites
    n.loc <- locus.params[1, 1]
    if(pop.data[1, 3] == "?") {
      full.seq <- paste(rep("A", n.loc), collapse = "")
      pop.data[, 3] <- rep(full.seq, nrow(pop.data))
    } else { # otherwise add A's to pad out to full sequence length
      partial.seq <- paste(rep("A", n.loc - nchar(pop.data[1, 3])), collapse = "")
      pop.data[, 3] <- sapply(pop.data[, 3], function(x) paste(x, partial.seq, sep = "", collapse = ""))
    }
    dna.seq <- strsplit(pop.data[, 3], "")
    names(dna.seq) <- pop.data[, 2]
    list(strata = data.frame(strata = pop.data[, 1]), dna.seq = as.DNAbin(dna.seq))
  } else {
    # compile diploid data
    n.loc <- ncol(pop.data) - 2
    pop.data <- do.call(rbind, lapply(seq(1, nrow(pop.data), 2), function(i) {
      ind <- pop.data[c(i, i + 1), ]
      locus.data <- as.vector(ind[, -(1:2)])
      c(ind[1, 1], paste(ind[, 2], collapse = "/"), locus.data)
    }))
    locus.data <- pop.data[, -c(1:2)]
    collapsed.loci <- do.call(cbind, lapply(seq(2, ncol(locus.data), by = 2), function(i) {
      a1 <- locus.data[, i - 1]
      a2 <- locus.data[, i]
      paste(a1, a2, sep = "/")
    }))
    colnames(collapsed.loci) <- paste("Locus", 1:ncol(collapsed.loci), sep = ".")
    rownames(collapsed.loci) <- pop.data[, 2]
    df2genind(collapsed.loci, sep = "/", pop = pop.data[, 1], type = "codom")
  }
  result
}


# run one iteration of fastsimcoal
fastsimcoal.skeleSim.run <- function(num.pops, Ne, sample.size = NULL, sample.time = NULL,
                                     growth.rate = NULL, mig.rates = NULL,
                                     hist.ev = NULL, num.chrom = 1,
                                     locus.params = NULL, label = "fastsimcoal.skeleSim",
                                     inf.site.model = TRUE, quiet = TRUE) {

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
  cmd <- paste("fastsimcoal -i", file, "-n 1",
               ifelse(inf.site.model, "-I", ""), ifelse(quiet, "-q", "")
  )
  err <- system(cmd, intern = F)

  if(err == 0) {
    if(!quiet) cat("fastsimcoal exited normally\n")
  } else {
    stop("fastsimcoal exited with error ", err, "\n")
  }

  invisible(file)
}


# wrapper for running one iteration of fastsimcoal and
#   reading results. return 'params' object has genetic data in
#   params$rep.result (dna = list of 2 elements, msat/snp = genind object)
sim.wrap.fastsimcoal <- function(params) {
  # create label for this run with scenario and replicate numbers
  current.label <- paste(
    params$label,
    params$common_params$current_scenario,
    params$common_params$current_replicate,
    sep = "."
  )
  params$common_params$current.label <- current.label

  # modify Ne and sample size if not sequence data to account for
  #   haploid -> diploid conversion of output
  size.mult <- if(params$common_params$locus_type == "sequence") 1 else 2

  # run one replicate of fastsimcoal
  fastsimcoal.skeleSim.run(
    num.pops = params$common_params$num_pops,
    Ne = params$common_params$pop_sizes * size.mult,
    sample.size = params$common_params$sample_sizes * size.mult,
    mig.rates = list(params$common_params$mig_rates),
    num.chrom = params$common_params$num_loci,
    hist.ev = params$fastsimcoal.params$hist.ev,
    sample.time = params$fastsimcoal.params$sample.time,
    growth.rate = params$fastsimcoal.params$growth.rate,
    locus.params = params$fastsimcoal.params$locus.params,
    inf.site.model = params$fastsimcoal.params$inf.site.model,
    label = current.label,
    quiet = params$quiet
  )

  arp.file <- paste(current.label, "_1_1.arp", sep = "")
  params$fastsimcoal.params$arp.file <- file.path(current.label, arp.file)
  params$rep.result <- fastsimcoal.skeleSim.read(params)
  params
}


# set parameters for fastsimcoal
set.fastsimcoal.params <- function(params) {
  params$fastsimcoal.params <- list()
  params$fastsimcoal.params$sample.times <- NULL
  params$fastsimcoal.params$growth.rate <- NULL
  params$fastsimcoal.params$inf.site.model <- TRUE

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
  locus.type <- switch(params$common_params$locus_type,
                       microsat = "MICROSAT", snp = "SNP", sequence = "DNA"
  )
  if(locus.type == "DNA") params$common_params$num_loci <- 1
  num.loci <- params$common_params$num_loci

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
