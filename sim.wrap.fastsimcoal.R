sim.wrap.fastsimcoal <- function(params) {
  params$fastsimcoal.params$last.file <- fastsimcoal.skeleSim.run(
    num.pops = params$common_params$num_pops,
    Ne = params$common_params$pop_sizes,
    sample.size = params$common_params$sample_sizes,
    mig.rates = list(params$common_params$mig_rates),
    num.chrom = params$common_params$num_loci,
    hist.ev = params$fastsimcoal.params$hist.ev,
    sample.time = params$fastsimcoal.params$sample.time,
    growth.rate = params$fastsimcoal.params$growth.rate,
    locus.params = params$fastsimcoal.params$locus.params,
    inf.site.model = params$fastsimcoal.params$inf.site.model,
    label = params$label, quiet = params$quiet
  )

  #invisible(fastsimcoal.skeleSim.read(params))
  params
}