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
fsc.run <- function(params) {
  label <- currentLabel(params)
  sc <- currentScenario(params)

  # Check that folder is empty
  if(file.exists(params@wd)) for(f in dir(label, full.names = T)) file.remove(f)

  locus.params <- sc@simulator.params@locus.params
  ploidy <- attr(locus.params, "ploidy")
  if(is.null(ploidy)) ploidy <- 1

  # Write fastsimcoal input file
  file <- fsc.write(
    num.pops = sc@num.pops,
    Ne = sc@pop.size,
    sample.size = sc@sample.size,
    sample.times = sc@simulator.params@sample.times,
    growth.rate = sc@simulator.params@growth.rate,
    mig.rates = sc@migration,
    hist.ev = sc@simulator.params@hist.ev,
    num.chrom = sc@simulator.params@num.chrom,
    locus.params = locus.params,
    ploidy = ploidy,
    label = label
  )

  # Run fastsimcoal
  num.cores <- params@num.cores
  cores.spec <- if(!is.null(num.cores)) {
    num.cores <- max(1, num.cores)
    num.cores <- min(num.cores, min(detectCores(), 12))
    paste(c("-c", "-B"), num.cores, collapse = " ")
  } else "-c 0 -B 12"
  cores.spec <- ""
  cmd.line <- paste(
    sc@simulator.params@fastsimcoal.exec, "-i", file, "-n 1",
    ifelse(params@quiet, "-q", ""), "-S", cores.spec
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

  num.chrom <- sc@simulator.params@num.chrom
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
fsc.write <- function(num.pops, Ne, sample.size = NULL, sample.times = NULL,
                      growth.rate = NULL, mig.rates = NULL,
                      hist.ev = NULL, num.chrom = NULL,
                      locus.params = NULL, ploidy = NULL, label = NULL) {

  opt <- options(scipen = 999)

  if(is.null(ploidy)) {
    pl <- attr(locus.params, "ploidy")
    ploidy <- if(is.null(pl)) 1 else pl
  }
  Ne <- Ne * ploidy
  if(!is.null(sample.size)) sample.size <- sample.size * ploidy

  if(is.null(label)) label <- "fastsimcoal.skeleSim"
  file <- paste(label, ".par", sep = "")
  mig.rates <- if(!is.null(mig.rates)) if(is.list(mig.rates)) mig.rates else list(mig.rates)
  hist.ev <- if(is.list(hist.ev)) do.call(rbind, hist.ev) else rbind(hist.ev)

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
fsc.read <- function(file, chrom.pos, ploidy) {
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

  # get data type
  data.type <- f[grep("DataType=", f)]
  data.type <- gsub("\tDataType=", "", data.type)
  switch(data.type,
    DNA = { # diploid SNPs
      dna.seq <- do.call(rbind, strsplit(data.mat[, 3], ""))
      if(all(dna.seq %in% 0:1)) {
        gen.data <- formatGenotypes(cbind(data.mat[, 1:2], dna.seq), ploidy)
        df2gtypes(gen.data, ploidy, description = file)
      } else { # haploid DNA sequences
        rownames(dna.seq) <- data.mat[, 2]
        dna.seq <- tolower(dna.seq)
        # create multidna object
        dna.seq <- new("multidna", lapply(1:nrow(chrom.pos), function(i) {
          as.matrix(dna.seq)[, chrom.pos[i, "start"]:chrom.pos[i, "end"]]
        }))
        # create gtypes object
        g <- sequence2gtypes(dna.seq, strata = data.mat[, 1], description = file)
        labelHaplotypes(g)$gtype
      }
    },
    MICROSAT = {
      gen.data <- formatGenotypes(data.mat, ploidy)
      # create gtypes object
      df2gtypes(gen.data, ploidy, description = file)
    }
  )
}

