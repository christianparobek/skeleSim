# wrapper for running one iteration of fastsimcoal and
#   reading results. return 'params' object has genetic data in
#   params$rep.result (dna = list of 2 elements, msat/snp = genind object)
fsc.run <- function(params) {
  # create label for this run with scenario and replicate numbers
  lbl <- currentLabel(params)

  # modify Ne and sample size if not sequence data to account for
  #   haploid -> diploid conversion of output
  size.mult <- if(params$common_params$locus_type == "sequence") 1 else 2

  # run one replicate of fastsimcoal
  fsc.read(
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
  arp.file <- paste(current.label, "_1_1.arp", sep = "")
  params$fastsimcoal.params$arp.file <- file.path(current.label, arp.file)
  params$rep.result <- fastsimcoal.skeleSim.read(params)
  params
}


# check that parameters are ready to be run
fsc.paramCheck <- function(params) {

  return(params)
}
