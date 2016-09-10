####### main shinyServer
####### sources most elements, though there are some
####### components at the bottom of this function
#######
shinyServer(function(input, output, session) {
#################### load some initializations that need to be within shinyServer
    source("serverInit.R", local = TRUE)
  ##################### parameter loading and saving
    source("saveParams.R", local = TRUE)
    source("loadParams.R", local = TRUE)

#    getssClass.file <- callModule(fileIO,"files",
#                                  roots=VolumeRoots,
#                                  objLabel = supportValues$objLabel,
#                                  ssLoadEnv = supportValues$ssLoadEnv,
#                                  ssClass = rValues$ssClass
#                                  )
#    rValues$ssClass <- getssClass.file()
  ##################### include the server code for Christians implemntation of
  ##################### the initial skelesim questions
    source("intro-questions-server.R",local=T)

  #################### general parameters for the simulation
  ####################

#    source("genparam-server.R",local=T)
    
  ##################### scenario helpers
  #################### stored in scenarios.R
    source("scenarios-server.R",local=T)

##################### simcoal helpers
#################### stored in simcoal-server.R
########################################
    source("simcoal-server.R",local=T)


################## rmetasim helpers
    source("rmetasim-server.R",local=T)
    
  ############## plotting
    source("serverplots.R",local=T)


  ##################### parameter loading and saving
    source("runSim.R", local = T)

######################## skeleSim class setup
    source("make-skelesim-class.R",local=T)

##########updating user interface
    source("update-ui.R",local=T)
    
##################visualization
source("vizServer.R", local = T)

#################sanity checks
source("sanity-checks.R", local=T)

####navbar header
    output$simtext <- renderText({
        if (is.null(rValues$ssClass@simulator.type))
            {"No Simulator"}
        else
            {
                ifelse(rValues$ssClass@simulator.type=="c","FastSimCoal params","Rmetasim Params")
            }
    })
    
####error messages.  Drop a message into rValues$msg to display to screen
    output$msg <- renderText({if (!is.null(rValues$msg)) rValues$msg})
#### quit button
    observeEvent(input$quitbtn,{stopApp()})
    
#############debugging
source("debug.R",local=T)

####### this needs to be replaced with Eric's loading code
#
#    observeEvent(input$readss,{
#        load("test.ssClass.rdata")
#        rValues$ssClass <- ssClass
#        rm(ssClass)
#    })


    
})

