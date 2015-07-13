# read a .arp file written by fastsimcoal
# to be called after simulation has been run but before analysis is to be run
fsc.read <- function(file, params) {
  stopifnot(require(adegenet))
  stopifnot(require(ape))
  f <- readLines(file)

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

  # return genetic data
  if(is.seq) { # if sequence data, return list with stratification and DNAbin sequences
    # replace sequence with all A's if there are no variable sites
    seq.len <- currentScenario(params)@sequence.length
    if(pop.data[1, 3] == "?") {
      full.seq <- paste(rep("A", seq.len), collapse = "")
      pop.data[, 3] <- rep(full.seq, nrow(pop.data))
    } else { # otherwise add A's to pad out to full sequence length
      partial.seq <- paste(rep("A", seq.len - nchar(pop.data[1, 3])), collapse = "")
      pop.data[, 3] <- sapply(pop.data[, 3], function(x) paste(x, partial.seq, sep = "", collapse = ""))
    }
    # format sequences to be converted to DNAbin ### !!!! CHANGE to APEX !!!!
    dna.seq <- lapply(strsplit(pop.data[, 3], ""), tolower)
    names(dna.seq) <- pop.data[, 2]
    list(strata = pop.data[, 1], dna.seq = as.DNAbin(dna.seq))
  } else { # if diploid data, return genind object
    n.loc <- ncol(pop.data) - 2
    # get population data
    pop.data <- do.call(rbind, lapply(seq(1, nrow(pop.data), 2), function(i) {
      ind <- pop.data[c(i, i + 1), ]
      locus.data <- as.vector(ind[, -(1:2)])
      c(ind[1, 1], paste(ind[, 2], collapse = "/"), locus.data)
    }))
    # extract allelic data
    locus.data <- pop.data[, -c(1:2)]
    # collapse loci into one column
    collapsed.loci <- do.call(cbind, lapply(seq(2, ncol(locus.data), by = 2), function(i) {
      a1 <- locus.data[, i - 1]
      a2 <- locus.data[, i]
      paste(a1, a2, sep = "/")
    }))
    colnames(collapsed.loci) <- paste("Locus", 1:ncol(collapsed.loci), sep = ".")
    rownames(collapsed.loci) <- pop.data[, 2]
    # create genind object
    df2genind(collapsed.loci, sep = "/", pop = pop.data[, 1], type = "codom")
  }
}
