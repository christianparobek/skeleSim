source("setup.R")

## !!!! CHECK IF coalParams CAN BE DELETED !!!!
#coalParams <<- new

ssClassInit <- function(){ #Just creates a skelesim class instance with one scenario
    ssClass <- new("skeleSim.params")
    ssClass@scenarios <- list(new("scenario.params"))
    ssClass@scenarios[[1]]@simulator.params <- new("fastsimcoal.params") #this will be changed if necessary reactively 
    ssClass
}

shinyServer(function(input, output,session) {

    rValues <- reactiveValues(ssClass=ssClassInit(),
                              scenarioNumber=1,
                              lstScenario=1,
                              history=NULL)
  ##################### parameter loading and saving
#  source("saveParams.R", local = TRUE)

  ##################### include the server code for Christians implemntation of
  ##################### the initial skelesim questions
  source("intro-questions-server.R",local=T)

  #################### general parameters for the simulation
  ####################

  source("genparam-server.R",local=T)
    
  ##################### scenario helpers
  #################### stored in scenarios.R
  source("scenarios-server.R",local=T)

  ##################### simcoal helpers
  #################### stored in simcoal-server.R
  ########################################
  source("simcoal-server.R",local=T)

  ############## plotting
  source("serverplots.R",local=T)


  ##################### parameter loading and saving
#  source("runSim.R", local = TRUE)

######################## skeleSim class setup
source("make-skelesim-class.R",local=T)


    
#############debugging
    output$ssClass <- renderTable({data.frame(item=(capture.output(str(rValues$ssClass))))})
####### this needs to be replaced with Eric's loading code
observeEvent(input$readss,{load("test.ssClass.rdata"); rValues$ssClass <- ssClass; rm(ssClass)})

    
})

