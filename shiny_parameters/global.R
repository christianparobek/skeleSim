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
source("../fastsimcoal/fsc.run.R")
source("../rmetasim/rms.run.R")

ssUserEnv <- new.env() # environment for holding session-persistent user objects
histry <- NULL     #saves a simcoal history
lstclick <- NULL    #last click
lstdblclick <- NULL #last double click
# coalParams <- new
ssClass <<- new("skeleSim.params")
objLabel <- NULL # syntactically valid name from 'title' slot of parameter object
fnameLabel <- NULL # combination of objLabel and timestamp for labelling filenames

