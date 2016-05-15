###
### implements a lot of the backend for the viz tab
###

###  selections set up when object read in

output$globalMainPlot <- renderPlot({

    gg.global(rValues$ssClass,stats=input$gstatsel,scenario=input$vizScenario)

})

output$locusMainPlot <- renderPlot({

    gg.locus(rValues$ssClass,stats=input$lstatsel,scenario=input$vizScenario)

})

output$pairwiseMainPlot <- renderPlot({
    
    gg.pairwise(rValues$ssClass,stats=input$pstatsel,scenario=input$vizScenario)

})
