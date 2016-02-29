#' @name fsc.run
#' @title Run fastsimcoal
#' @description Run fastsimcoal
#'
#' @param params a \linkS4class{skeleSim.params} object.
#'
#' @return a modified \linkS4class{skeleSim.params} object with the results of
#'   a fastsimcoal run.
#'
#' @importFrom parallel detectCores
#' @export
#'
fsc.run <- function(pop.info, locus.params, mig.rates = NULL, hist.ev = NULL, label = NULL,
                    fsc.exec = "fsc252", quiet = TRUE, num.cores = NULL) {

  if(is.null(label)) label <- "fsc.run"
  if(file.exists(label)) for(f in dir(label, full.names = T)) file.remove(f)

  # Write fastsimcoal input file
  file <- fsc.write(
    pop.info = pop.info, locus.params = locus.params,
    mig.rates = mig.rates, hist.ev = hist.ev, label = label
  )

  # Run fastsimcoal
  cores.spec <- if(!is.null(num.cores)) {
    num.cores <- max(1, num.cores)
    num.cores <- min(num.cores, min(detectCores(), 12))
    paste(c("-c", "-B"), num.cores, collapse = " ")
  } else ""
  cmd.line <- paste(
    fsc.exec, "-i", file, "-n 1",
    ifelse(quiet, "-q", ""), "-S", cores.spec,
    attr(locus.params, "opts")
  )
  err <- if(.Platform$OS.type == "unix") {
    system(cmd.line, intern = F)
  } else {
    shell(cmd.line, intern = F)
  }

  if(err == 0) {
    if(!quiet) cat("fastsimcoal exited normally\n")
  } else {
    stop("fastsimcoal exited with error ", err, "\n")
  }

  arp.file <- file.path(label, paste(label, "_1_1.arp", sep = ""))
  fsc.read(arp.file, locus.params)
}


#' @rdname fsc.run
#'
#' @param num.pops number of populations.
#' @param Ne numeric vector: size of each population.
#' @param sample.size numeric vector: number of samples to draw from each population.
#' @param sample.times numeric vector: time step samples are drawn from for each population.
#' @param growth.rate numeric vector: growth rate of each population.
#' @param mig.rates list of matrices: migration rates between populations.
#' @param hist.ev matrix: one row per historical event.
#' @param num.chrom number of unique chromosomes.
#' @param locus.params data.frame: one row per locus parameter.
#' @param label character: label for file naming.
#'
fsc.write <- function(pop.info, locus.params, mig.rates = NULL, hist.ev = NULL, label = NULL) {

  opt <- options(scipen = 999)

  ploidy <- attr(locus.params, "ploidy")
  pop.info[, c("pop.size", "sample.size")] <- pop.info[, c("pop.size", "sample.size")] * ploidy

  if(is.null(label)) label <- "fastsimcoal.output"
  file <- paste(label, ".par", sep = "")
  mig.rates <- if(!is.null(mig.rates)) if(is.list(mig.rates)) mig.rates else list(mig.rates)
  hist.ev <- if(is.list(hist.ev)) do.call(rbind, hist.ev) else rbind(hist.ev)

  # Write input file
  write(paste("//  <<", label, ">>  (input from 'fastsimcoal.skeleSim.run')"), file)
  write(paste(nrow(pop.info), "populations to sample"), file, append = T)

  write("//Population effective sizes", file, append = T)
  for(i in 1:nrow(pop.info)) write(pop.info[i, "pop.size"], file, append = T)

  write("//Sample sizes", file, append = T)
  for(i in 1:nrow(pop.info)) write(pop.info[i, c("sample.size", "sample.times")], file, append = T)

  write("//Growth rates", file, append = T)
  for(i in 1:nrow(pop.info)) write(pop.info[i, "growth.rate"], file, append = T)

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

  num.chrom <- attr(locus.params, "num.chrom")
  if(!is.null(num.chrom)) locus.params$chromosome <- 1
  locus.params <- split(locus.params, locus.params$chromosome)
  num.independent <- if(is.null(num.chrom)) length(locus.params) else num.chrom
  chrom.struct <- if(num.independent == 1 | !is.null(num.chrom)) 0 else 1
  write("//Number of independent loci [chromosomes]", file, append = T)
  write(paste(num.independent, chrom.struct), file, append = T)
  for(block in locus.params) {
    block$chromosome <- NULL
    write("//Per chromosome: Number of linkage blocks", file, append = T)
    write(nrow(block), file, append = T)
    write("//Per block: data type, num loci, rec. rate and mut rate + optional parameters", file, append = T)
    for(i in 1:nrow(block)) {
      line <- paste(block[i, ], collapse = " ")
      write(gsub(" NA", "", line), file, append = T)
    }
  }

  options(opt)
  invisible(file)
}


#' @rdname fsc.run
#'
#' @param file character filename
#' @param chrom.pos a matrix giving the start and end positions of each chromosome.
#' @param ploidy a numeric giving the ploidy of the loci (1 = haploid,
#'   2 = diploid).
#'
#' @import strataG
#' @importFrom stringi stri_extract_last_regex
#' @importFrom swfscMisc zero.pad
#'
fsc.read <- function(file, locus.params) {
  formatGenotypes <- function(x, ploidy) {
    # reformat matrix to have alleles side-by-side
    nloci <- ncol(x) - 2
    loc.end <- seq(ploidy, nrow(x), by = ploidy)
    gen.data <- do.call(rbind, lapply(loc.end, function(i) {
      allele.i <- (i - ploidy + 1):i
      loci <- as.vector(x[allele.i, -(1:2)])
      id <- paste(x[allele.i, 2], collapse = ".")
      pop <- x[allele.i[1], 1]
      c(id, pop, loci)
    }))
    # rename loci
    locus_names <- paste("Locus", zero.pad(1:nloci), sep = "_")
    locus_names <- paste(rep(locus_names, each = ploidy), 1:ploidy, sep = ".")
    colnames(gen.data) <- c("id", "pop", locus_names)
    gen.data
  }

  formatDNA <- function(dna.seq, pop, locus.params) {
    # create multidna object splitting chromosomes into loci
    num.chrom <- attr(locus.params, "num.chrom")
    chrom.pos <- if(is.null(num.chrom)) {
      tapply(locus.params$num.markers, locus.params$chromosome, sum)
    } else {
      rep(sum(locus.params$num.markers), num.chrom)
    }
    chrom.pos <- cumsum(chrom.pos)
    chrom.pos <- cbind(start = c(1, chrom.pos[-length(chrom.pos)] + 1), end = chrom.pos)

    rownames(dna.seq) <- pop
    dna.seq <- tolower(dna.seq)
    new("multidna", lapply(1:nrow(chrom.pos), function(i) {
      as.matrix(dna.seq)[, chrom.pos[i, "start"]:chrom.pos[i, "end"]]
    }))
  }

  f <- readLines(file)

  # get start and end points of data blocks
  start <- grep("SampleData=", f) + 1
  end <- which(f == "}") - 2
  pos <- cbind(start, end)

  # compile data into 3 column character matrix
  data.mat <- do.call(rbind, lapply(1:nrow(pos), function(i) {
    f.line <- f[pos[i, 1]:pos[i, 2]]
    f.line <- gsub("[[:space:]]+", "--", f.line)
    result <- do.call(rbind, strsplit(f.line, "--"))[, -2]
    cbind(rep(paste("Sample", i), nrow(result)), result)
  }))

  ploidy <- attr(locus.params, "ploidy")

  # get data type
  data.type <- f[grep("DataType=", f)]
  data.type <- gsub("\tDataType=", "", data.type)
  switch(
    data.type,
    DNA = { # diploid SNPs
      dna.seq <- do.call(rbind, strsplit(data.mat[, 3], ""))
      if(attr(locus.params, "ploidy") == 2) {
        gen.data <- formatGenotypes(cbind(data.mat[, 1:2], dna.seq), ploidy)
        df2gtypes(gen.data, ploidy, description = file)
      } else { # haploid DNA sequences
        dna.seq <- formatDNA(dna.seq, data.mat[, 2], locus.params)
        g <- sequence2gtypes(dna.seq, strata = data.mat[, 1], description = file)
        labelHaplotypes(g)$gtype
      }
    },
    MICROSAT = {
      gen.data <- formatGenotypes(data.mat, ploidy)
      df2gtypes(gen.data, ploidy, description = file)
    },
    NULL
  )
}


fsc.popInfo <- function(pop.size, sample.size, sample.times = 0, growth.rate = 0) {
  cbind(
    pop.size = pop.size, sample.size = sample.size,
    sample.times = sample.times, growth.rate = growth.rate
  )
}

fsc.locusParams <- function(locus.type = c("dna", "msat", "snp"),
  sequence.length = NULL, num.loci = NULL, mut.rate = NULL,
  transition.rate = 1 / 3, gsm.param = 0,
  range.constraint = 0, recomb.rate = 0, chromosome = NULL,
  num.chrom = NULL, ploidy = NULL) {

  createLocusParams <- function(chr, type, num.markers, recomb.rate, param.4,
                                param.5, param.6, ploidy, num.chrom) {
    if(is.null(chr)) chr <- 1
    df <- data.frame(
      chromosome = chr, type = type, num.markers = num.markers,
      recomb.rate = recomb.rate, param.4 = param.4, param.5 = param.5,
      param.6 = param.6, stringsAsFactors = FALSE
    )
    df <- df[order(df$chromosome), ]
    attr(df, "num.chrom") <- num.chrom
    attr(df, "ploidy") <- ploidy
    return(df)
  }

  df <- switch(
    match.arg(locus.type),
    dna = createLocusParams(
      chromosome, "DNA", sequence.length, recomb.rate, mut.rate,
      transition.rate, NA, 1, num.chrom
    ),
    msat = createLocusParams(
      chromosome, "MICROSAT", num.loci, recomb.rate, mut.rate, gsm.param,
      range.constraint, 2, num.chrom
    ),
    snp = createLocusParams(
      chromosome, "DNA", 1, recomb.rate, mut.rate, 1,
      NA, 2, num.loci
    )
  )
  attr(df, "opts") <- if(locus.type == "snp") "-s" else ""
  return(df)
}


#' @title Load skeleSim scenario parameters for fastsimcoal
#' @description Load skeleSim scenario parameters for fastsimcoal
#'
#' @param num.pops number of populations.
#' @param pop.size a vector giving size of each populaiton.
#' @param sample.size a vector giving the number of samples to take from each
#'   population.
#' @param migration a \code{num.pop} x \code{num.pop} matrix or list of matrices
#'   giving the migration rates between each population.
#' @param locus.type a character representation of what type of marker to simulate.
#'   Can be "dna", "msat", or "snp".
#' @param num.loci \code{msat, snp}: number of loci to simulate.
#' @param sequence.length \code{dna}: number of DNA base pairs to use.
#' @param mut.rate \code{dna, msat}: per base pair or locus mutation rate.
#' @param sample.times a vector giving the number of generations in the past
#'   at which samples are taken.
#' @param growth.rate a vector giving the growth rate of each population.
#' @param hist.ev a matrix describing historical events.
#' @param num.chrom a value giving the number of chromosomes that the
#'   \code{locus.params} marker specifications should be copied for. If
#'   \code{NULL}, then chromosome assignment is taken from the
#'   \code{chromosome} column in \code{locus.params}. Any non-\code{NULL}
#'   integer will cause the \code{chromosome} column to be ignored.
#' @param transition.rate dna: fraction of substitutions that are transitions.
#' @param recomb.rate recombination rate between adjacent markers.
#' @param chromosome number or character identifying which chromosome the marker
#'   is on.
#' @param min.freq \code{snp}: minimum frequency for the derived allele.
#' @param gsm.param \code{msat}: Value of the geometric parameter for a
#'   Generalized Stepwise Mutation (GSM) model. This value represents the
#'   proportion of mutations that will change the allele size by more than
#'   one step. Values between 0 and 1 are required. A value of 0 is for a
#'   strict Stepwise Mutation Model (SMM).
#' @param range.constraint \code{msat}: Range constraint (number of different
#'   alleles allowed). A value of 0 means no range constraint.
#'
#' @note Vectors for \code{pop.size, sample.size, sample.times, and growth.rate}
#'   will be expanded/recycled to ensure they are as long as \code{num.pops}.
#'
#' Depending on the choice of \code{locus.type}, values for some arguments may
#'   be ignored. See argument list above for which arguments are applicable
#'   to which \code{locus.type}.
#'
#'
#'
#' @return a \linkS4class{scenario.params} object to be loaded into a list in the
#'   \code{scenarios} slot of a \linkS4class{skeleSim.params} object.
#'
#' @export
#'


#' @title Create fastsimcoal historical event matrices
#' @description Create fastsimcoal historical event matrices
#'
#' @param num.gen Number of generations, t, before present at which the
#'   historical event happened.
#' @param source.deme Source deme (the first listed deme has index 0)
#' @param sink.deme Sink deme
#' @param prop.migrants Expected proportion of migrants to move from source to sink.
#' @param new.sink.size New size for the sink deme, relative to its size at
#'   generation t.
#' @param new.sink.growth New growth rate for the sink deme.
#' @param new.mig.mat New migration matrix to be used further back in time.
#' @param hist.ev a matrix describing historical events, with one row per event.
#' @param pop.size numerical vector giving size of each population.
#' @param growth.rate numerical vector giving growth rate of each population.
#' @param num.mig.mats number of migration matrices.
#'
#' @return a blank fastsimcoal historical event matrices that can be
#'   filled in later
#'
#' @export
#'
fsc.histEvMat <- function(num.gen = 0, source.deme = 0, sink.deme = 0,
                          prop.migrants = 1, new.sink.size = 1,
                          new.sink.growth = 0, new.mig.mat = 0) {
  cbind(
    num.gen = num.gen, source.deme = source.deme, sink.deme = sink.deme,
    prop.migrants = prop.migrants, new.sink.size = new.sink.size,
    new.sink.growth = new.sink.growth, new.mig.mat = new.mig.mat
  )
}