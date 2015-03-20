############### Take test data FASTA file and trial ########################



DNAbin2gtypes <- function(go, strata, id.pop = NULL){

  if(length(stratafoo > 1)){
    sv <- which(names(strata) == id.pop)
  } else {
    sv <- 1
  }

  go <- as.matrix(go)

  x <- sapply(rownames(go), function(n) {
    as.character(go[n, ])[1, ]
  }, simplify = FALSE)

  gtypes(strata, strata.col = sv, dna.seq = x)

}
