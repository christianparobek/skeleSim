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

histry <- NULL     #saves a simcoal history
lstclick <- NULL    #last click
lstdblclick <- NULL #last double click
# coalParams <- new
#ssClass <<- new("skeleSim.params")
fnameLabel <- NULL # combination of ssClass@label and timestamp for labelling filenames