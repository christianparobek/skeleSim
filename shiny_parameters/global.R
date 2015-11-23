# load some packages
# source some files
#
options(shiny.trace = F)  # change to T for trace
require(shiny)
require(shinyFiles)
require(shinyIncubator)
require(igraph)

source("utils.R")
source("mig.matrix.R")
source("simcoal-history.R")
source("../skeleSim.classes.R")
source("../fastsimcoal/fsc.run.R")
source("../fastsimcoal/fsc.classes.R")
source("../fastsimcoal/fsc.scenarioCheck.R")
source("../rmetasim/rms.run.R")
source("../rmetasim/rms.classes.R")
source("../rmetasim/rms.scenarioCheck.R")


###these two globals record a record of the last
###click on the history graphs
lstclick <- NULL    #last click
lstdblclick <- NULL #last double click





