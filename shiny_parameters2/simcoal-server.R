#
# file to be included in server.R that specifies parts of the simcoal backend
#

hst <- reactive({
    if (is.null(histry))
        {
            histry<<-create.new.history(npop=input$numpops)
        }  else  {
            histry<<-simcoal.history.change(histry,input$histplotClick,
                                             input$histplotDblclick)
            print(input$histplotClick)
        }
    return(histry)
})

output$simhistPlot <- renderPlot({
    simcoal.history.plot(hst())
})

output$simhistTbl <- renderTable({
    hst()
})
