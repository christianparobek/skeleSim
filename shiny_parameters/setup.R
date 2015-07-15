#
# load some packages
# source some files
#
options(shiny.trace = F)  # cahnge to T for trace
require(shiny)
require(igraph)

source("utils.R")
source("mig.matrix.R")
source("simcoal-history.R")
source("../skeleSim.classes.R")
source("../skeleSim.funcs.R")
source("../rmetasim/rms.run.R")
source("../fastsimcoal/fsc.run.R")
