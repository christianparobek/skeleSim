# AES 7/2/15
#
# these are renderUI statements that only make sense
# included into the skelesim parameter interface
# using source()

#scenario UIs  this allows the defaults to be set for particular scenarios

output$scenarioNumberUI <- renderUI({
    numericInput("scenarioNumber", "Which scenario",
             value = 1)
})

output$mutrateUI <- renderUI({
    numericInput("mut.rate", "Mutation Rate",
             value = 1e-4)
})



output$mutrateUI <- renderUI({
    numericInput("mut.rate", "Mutation Rate",
             value = 1e-4)
})


output$migmodelUI <- renderUI({
selectInput("migModel", "Migration Model",
            c("island","stepping.stone.linear",
              "stepping.stone.circular","twoD","twoDwDiagonal","distance"))
})

output$migrateUI <- renderUI(
    {numericInput("migRate", "Migration rate",1)}
    )
                                        #
output$rows <- renderUI({
    if (input$migModel%in%c("twoD","twoDwDiagonal","distance"))
        {
            numericInput("rows", "Rows in a grid-shaped landscape",2)
        }
})

output$cols <- renderUI({
    if (input$migModel%in%c("twoD","twoDwDiagonal","distance"))
        {
            numericInput("cols", "Cols in a grid-shaped landscape",2)
        }
})


output$distanceFun <- renderUI({
    if (input$migModel%in%c("distance"))
        {
            textInput("distfun", "Distance function (must be an R function)","dexp")
        }
})


################SIMCOAL SPECIFICS
output$infsitesUI <- renderUI({
    if (input$coalescent==T)
        {
            checkboxInput("infSiteModel", "Infinite site model",
                          value = FALSE)
        } else {
        }
})        

output$simhistUI <- renderUI({
    if (input$coalescent==T)
        {
            checkboxInput("simhistflag", "Specify a demographic history for simulation?",
                          value = FALSE)
        } else {
        }
})        


