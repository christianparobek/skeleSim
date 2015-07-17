source("setup.R")

## !!!! CHECK IF coalParams CAN BE DELETED !!!!
#coalParams <<- new

shinyServer(function(input, output,session) {
  ##################### parameter loading
  source("paramsLoad.R", local = TRUE)

  ##################### parameter saving
  source("paramsSave.R", local = TRUE)

  ##################### running simulator
  source("simRun.R", local = TRUE)

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
  source("make-skelesim-class.R",local=T)

  #############debugging
  output$tbl <- renderTable({
    inmat()
  })
  output$txt <- renderText({
    #    names(input)
    #    c(input$coalescent,input[["coalescent"]])
  })

})
