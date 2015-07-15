source("setup.R")

####
### dreaded global variables
histry <<- NULL     #saves a simcoal history
lstclick <<- NULL    #last click
lstdblclick <<- NULL #last double click
scenarios <<- vector("list",1)   #will be a list of scenarios, start with one


shinyServer(function(input, output,session) {

##################### include the server code for Christians implemntation of
##################### the initial skelesim questions
source("intro-questions-server.R",local=T)

##################### server-side user interface specifcations
#####################  are in the file renderUI.R
source("renderUI.R",local=T)

##################### scenario helpers
#################### stored in scenarios.R
source("scenarios.R",local=T)

##################### simcoal helpers
#################### stored in simcoal-server.R
########################################
source("simcoal-server.R",local=T)

############## plotting
source("serverplots.R",local=T)

######################## skeleSim class setup
#source("make-skelesim-class.R",local=T)

#############debugging
                output$tbl <- renderTable({
                    inmat()
                })
                output$txt <- renderText({
                                        #    names(input)
                                        #    c(input$coalescent,input[["coalescent"]])
                })

  observeEvent(input$btnRunSim, {
    # check parameters

    # write files
    label <- "test"
    fname <- gsub("[[:punct:]]", ".", label)
    paramsFname <- paste(fname, ".rdata", sep = "")
    save(label, file = paramsFname)
    scriptFname <- paste(fname, ".R", sep = "")
    write("rm(list = ls())", file = scriptFname)
    write("load('test.params.ws.rdata')", file = scriptFname, append = TRUE)
    write("test.params <- runSim(test.params)", file = scriptFname, append = TRUE)
    outFname <- paste(fname, "Output.rdata", sep = "")
    line <- paste("save(test.params, file = '", outFname, "')", sep = "")
    write(line, file = scriptFname, append = TRUE)
    output$runText <- renderText({
      paste("Script file: '", scriptFname, "' written.", sep = "")
    })
    # run system command to execute script file

  })
})
