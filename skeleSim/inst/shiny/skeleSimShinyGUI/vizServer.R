###
### implements a lot of the backend for the viz tab
###

###  selections set up when object read in

output$globalMainPlot <- renderPlot({

    gg.global(rValues$ssClass,stats=input$gstatsel)

})

output$locusMainPlot <- renderPlot({

    gg.locus(rValues$ssClass,stats=input$lstatsel)

})

output$pairwiseMainPlot <- renderPlot({

    gg.pairwise(rValues$ssClass,stats=input$pstatsel)

})
