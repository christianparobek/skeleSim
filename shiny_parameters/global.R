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

####this global list contains values that are needed for file operations
####cannot use reactive due to constraints imposed by shinyFiles
supportValues <- list(ssLoadEnv=new.env(),  #environment to load an rdata file into
                      objLabel=NULL,        #name of ssClass object for saving
                      roots = c(getVolumes()(),home="~",temp=tempdir(),wd="./"),   #function from shinyFiles
                      simroot = NULL
                      )




