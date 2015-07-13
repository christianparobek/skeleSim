# AES 7/2/15
#
# these are renderUI statements that only make sense
# included into the skelesim parameter interface
# using source()
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


