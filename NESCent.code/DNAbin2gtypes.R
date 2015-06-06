############### Take test data FASTA file and trial ########################
## takes a FASTA file


DNAbin2gtypes <- function(go, strata, id.pop = 1){

  sv <- id.pop


  go <- as.matrix(go)

  x <- sapply(rownames(go), function(n) {
    as.character(go[n, ])[1, ]
  }, simplify = FALSE)

  gtypes(strata, strata.col = sv, dna.seq = x)

}
