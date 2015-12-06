# load some packages
# source some files
#
#options(shiny.trace = F)  # change to T for trace
#require(shiny)
#require(shinyFiles)
#require(shinyIncubator)
#require(igraph)

source("utils.R")
source("mig.matrix.R")
source("simcoal-history.R")
#source("../skeleSim.classes.R")
#source("../fastsimcoal/fsc.run.R")
#source("../fastsimcoal/fsc.classes.R")
#source("../fastsimcoal/fsc.scenarioCheck.R")
#source("../rmetasim/rms.run.R")
#source("../rmetasim/rms.classes.R")
#source("../rmetasim/rms.scenarioCheck.R")


###these two globals record a record of the last
###click on the history graphs
lstclick <- NULL    #last click
lstdblclick <- NULL #last double click



####### helper functions
#######
ssClassInit <- function(){ #Just creates a skelesim class instance with one scenario
    ssClass <- new("skeleSim.params")

##### just a placeholder !!!!!!!! #######
##### need to include actual rep analyses
        ssClass@rep.analysis.func <-  function(params) {
            result = rnorm(5)
            names(result) <- paste("result", 1:5, sep = ".")
            params@rep.result <- result
            params
        }
##############################################
    
    ssClass@scenarios <- list(new("scenario.params"))

    #default values for ssClass@scenarios  (could be set in class definition)
    ssClass@scenarios[[1]]@num.pops <- 1
    ssClass@scenarios[[1]]@pop.size <- 100
    ssClass@scenarios[[1]]@sample.size <- 10
    ssClass@scenarios[[1]]@migration <- list(matrix(0,nrow=1,ncol=1))
    ssClass@scenarios[[1]]@mig.helper <-
        list(migModel="island",migRate=1,rows=1,cols=1,distfun="dexp")
    ssClass@scenarios[[1]]@num.loci <- 1
    ssClass@scenarios[[1]]@sequence.length <- 100
    ssClass@scenarios[[1]]@mut.rate <- 10e-5
    ssClass@scenarios[[1]]@simulator.params <-
        new("fastsimcoal.params") #this will be changed if necessary reactively

    ssClass
}




