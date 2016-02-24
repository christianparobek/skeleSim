#' @name fsc.read
#' @title Read fastsimcoal output file
#' @description Read fastsimcoal output file
#'
#' @param file character filename
#' @param chrom.pos a matrix giving the start and end positions of each chromosome.
#' @param ploidy a numeric giving the ploidy of the loci (1 = haploid,
#'   2 = diploid).
#'
#' @return genetic data from \code{file} as a \linkS4class{gtypes} object.
#'
#' @import strataG
#' @importFrom stringi stri_extract_last_regex
#'
fsc.read <- function(file, chrom.pos, ploidy) {
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
    DNA = {
      dna.seq <- do.call(rbind, strsplit(data.mat[, 3], ""))
      if(all(dna.seq %in% 0:1)) {
        gen.data <- formatGenotypes(cbind(data.mat[, 1:2], dna.seq), ploidy)
        df2gtypes(gen.data, ploidy, description = file)
      } else {
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

#' @rdname fsc.read
#'
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
  locus_names <- paste("Locus", 1:nloci, sep = "_")
  locus_names <- paste(rep(locus_names, each = ploidy), 1:ploidy, sep = ".")
  colnames(gen.data) <- c("id", "pop", locus_names)
  gen.data
}
