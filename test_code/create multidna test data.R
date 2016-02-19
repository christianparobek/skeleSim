library(strataG)
data("dolph.strata")
data("dolph.seqs")
dolph.seqs <- mafft(dolph.seqs)
fine <- dolph.strata$fine
names(fine) <- dolph.strata$id

dna.seqs <- as.multidna(lapply(list(1:150, 151:300, 301:402), function(i) {
  as.DNAbin(as.character(as.matrix(dolph.seqs))[, i])
}))

g <- labelHaplotypes(sequence2gtypes(dna.seqs, strata = fine))$gtypes
