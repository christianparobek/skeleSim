#' @title Create fastsimcoal locus parameter matrices
#' @description Create blank fastsimcoal locus parameter matrices that can be
#'   filled in later
#'
#' @param sc a \linkS4class{scenario.params} object.

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
