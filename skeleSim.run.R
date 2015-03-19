rm(list = ls())
gc()
source("fastsimcoal.skeleSim.run.r")
source("genind.metadata.getter.r")
source("set.commonparams.r")
source("set.fastsimcoal.params.r")
source("sim.choice.r")
source("sim.wrap.fastsimcoal.r")
source("fastsimcoal.skeleSim.read.r")
library(rmetasim)

params <- list(
  label = "fastsimcoal.test",
  quiet = FALSE,
  user_has_data = FALSE
)
params <- set.commonparams(params, NULL)


params <- set.fastsimcoal.params(params)

params$