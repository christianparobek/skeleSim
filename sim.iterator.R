## This is a function to iterate through simulations and analyses
## This function should take the following information:
##      The simulation parameters contained in the megalist
##      The analyis function contained in the megalist
## And this function should return:
##      Modifications to the megalist in the replicate metadata
##      Mods to the megalist, adding the output of the analyses

## Started 19 March 2015
## Started by cp

############################
############################

## Load required packages
library(adegenet)

## Load required functions (that we wrote)
source("new.mainparam.list.R")
source("set.commonparams.R")

## Make a toy megalist
mainparams_test <- new.mainparam.list()
set.commonparams(mainparams_test)

## Make a toy analysis function


#############################
#### Run the simulations ####
#############################

sim.iterator <- function(megalist){
  
  simfun <- 
  analysisfun <- 
}