source("setup.R")

####
### dreaded global variable
histry <<- NULL  #saves a simcoal history
lstclick <- NULL
lstdblclick <- NULL
shinyServer(function(input, output,session) {

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
                
            })
