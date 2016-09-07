# load some packages
# source some files
#
#options(shiny.trace = F)  # change to T for trace
#require(shiny)
#require(shinyFiles)
#require(shinyIncubator)
#require(igraph)
library(dplyr)
library(ggplot2)
source("utils.R")
source("mig.matrix.R")
source("simcoal-history.R")
source("vizAbstract.R")
source("vizfuncs.R")
#source("../skeleSim.classes.R")
#source("../fastsimcoal/fsc.run.R")
#source("../fastsimcoal/fsc.classes.R")
#source("../fastsimcoal/fsc.scenarioCheck.R")
#source("../rmetasim/rms.run.R")
#source("../rmetasim/rms.classes.R")
#source("../rmetasim/rms.scenarioCheck.R")
source("modules/matrixInputs.R")
source("modules/tableinput.R")


###these two globals record a record of the last
###click on the history graphs
#lstclick <- NULL    #last click
#lstdblclick <- NULL #last double click



####### helper functions
#######
ssClassInit <- function(){ #Just creates a skelesim class instance with one scenario
    ssClass <- new("skeleSim.params")

##### just a placeholder !!!!!!!! #######
##### need to include actual rep analyses
#        ssClass@rep.analysis.func <-  function(params) {
#            result = rnorm(5)
#            names(result) <- paste("result", 1:5, sep = ".")
#            params@rep.result <- result
#            params
#        }
##############################################
    ssClass@rep.analysis.func <- analysisFunc

#    ssClass@simulator.type <- "f"

    ssClass@scenarios <- list(new("scenario.params"))

    #default values for ssClass@scenarios  (could be set in class definition)
    ssClass@scenarios[[1]]@num.pops <- 2
    ssClass@scenarios[[1]]@pop.size <- c(100,100)
    ssClass@scenarios[[1]]@sample.size <- c(10,10)
    ssClass@scenarios[[1]]@migration <- list(matrix(c(0, 0.1, 0.1, 0),nrow=2,ncol=2))
    ssClass@scenarios[[1]]@mig.helper <-
        list(migModel="island",migRate=1,rows=1,cols=1,distfun="dexp")
    ssClass@scenarios[[1]]@num.loci <- 1
    ssClass@scenarios[[1]]@sequence.length <- 100
    ssClass@scenarios[[1]]@mut.rate <- 1e-5
    ssClass@current.scenario <- 1
    ssClass@current.replicate <- 1
    ssClass@scenarios[[1]]@simulator.params <-
        fastsimcoalInit(2) #2 is the number of demes. this will be changed if necessary reactively

    ssClass
}

fastsimcoalInit <- function(np){ #np num populations
    parms <- new("fastsimcoal.params")
    parms@growth.rate <- rep(0,np)
    parms@sample.times <- rep(0L,np) #must be integer
    parms
}

rmetasimInit <- function(np){ #np num pops
    parms <- new("rmetasim.params")
    parms@num.stgs <- 2
    parms@selfing <- 0
    parms@init.pop.sizes <- rep(c(100,100),np)
    parms@surv.matr <- matrix(c(0.2,0,
                               0.4,0.1),nrow=2,byrow=T)
    parms@repr.matr <- matrix(c( 0 ,4,
                               0 ,0),nrow=2,byrow=T)
    parms@male.matr <- matrix(c(0,0,
                               0,1),nrow=2,byrow=T)
    parms
}




