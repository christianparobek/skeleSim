#source("setup.R")

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

shinyServer(function(input, output,session) {
    
    VolumeRoots = c(home="~",getVolumes()(),temp=tempdir(),wd="./")   #function from shinyFiles

    rValues <- reactiveValues(ssClass=ssClassInit(),
                              scenarioNumber=1,
                              lstScenario=1,
                              history=NULL,
                              msg=NULL)

####this reactiveValue contains values that are needed for file operations
####cannot use reactive due to constraints imposed by shinyFiles
supportValues <- reactiveValues(ssLoadEnv=new.env(),  #environment to load an rdata file into
                                objLabel=NULL,        #name of ssClass object for saving
                                simroot = NULL
                                )

  ##################### parameter loading and saving
    source("saveParams.R", local = TRUE)
    source("loadParams.R", local = TRUE)

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
    source("runSim.R", local = TRUE)

######################## skeleSim class setup
    source("make-skelesim-class.R",local=T)

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

    observeEvent(input$quitbtn,{stopApp()})
    
#############debugging
    output$ssClass <- renderTable({data.frame(item=(capture.output(str(rValues$ssClass))))})

####### this needs to be replaced with Eric's loading code
#
#    observeEvent(input$readss,{
#        load("test.ssClass.rdata")
#        rValues$ssClass <- ssClass
#        rm(ssClass)
#    })


    
})

