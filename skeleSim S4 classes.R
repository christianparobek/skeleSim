setClass("skeleSim.params", slots = c(
  title = "character",
  date = "POSIXlt",
  question = "character", # (n)ull, p(o)wer, (p)erformance
  simulator = "character", # (c)oalescent or (f)orward-time
  user.has.data = "logical",
  scenarios = "list", # a list of scenario.params S4 class objects
  param.check.func = "function",
  sim.func = "function",
  last.sample = "matrix",
  rep.analysis.func = "function",
  rep.analysis = "numeric",
  sim.summary.func = "function",
  analysis.results = "matrix"
))

setClass("scenario.params", slots = c(
  num.pops = "integer",
  pop.size = "integer",
  sample.size = "integer",
  migration = "list",
  simulator.params = "ANY" # can be list or class of simulator-specific parameters
))

setClass("fastsimcoal.params", slots = c(
  sample.times = "numeric",
  growth.rate = "numeric",
  hist.ev = "matrix",
  locus.params = "list"
))
