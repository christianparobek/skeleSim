fastsimcoal.skeleSim.sequence.read <- function(file) {
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
  is.haploid <- switch(data.type, DNA = T, MICROSAT = F, STANDARD = F)
  if(is.haploid) {
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
    g <- gtypes(dna.seq, strata = pop.data[, 1], description = file)
    label.haplotypes(g, "Hap.")
  } else {
    # compile diploid data
    n.loc <- ncol(pop.data) - 2
    pop.data <- do.call(rbind, lapply(seq(1, nrow(pop.data), 2), function(i) {
      ind <- pop.data[c(i, i + 1), ]
      locus.data <- as.vector(ind[, -(1:2)])
      c(ind[1, 1], paste(ind[, 2], collapse = "/"), locus.data)
    }))
    pop.data <- data.frame(pop.data[, 2], pop.data[, 1], pop.data[, -(1:2)])
    gtypes(pop.data, description = file)
  }
}

fastsimcoal.skeleSim.read <- function(params) {




