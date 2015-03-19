sim.wrap.fastsimcoal <- function(params) {
  # create label for this run with scenario and replicate numbers
  current.label <- paste(
    params$label,
    params$common_params$current_scenario,
    params$common_params$current_replicate,
    sep = "."
  )
  params$common_params$current.label <- current.label

  # modify Ne and sample size if not sequence data to account for
  #   haploid -> diploid conversion of output
  size.mult <- if(params$common_params$locus_type == "sequence") 1 else 2

  # run one replicate of fastsimcoal
  fastsimcoal.skeleSim.run(
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

  arp.file <- paste(current.label, "_1_1.arp", sep = "")
  params$fastsimcoal.params$arp.file <- file.path(current.label, arp.file)
  params$rep.result <- fastsimcoal.skeleSim.read(params)
  params
}