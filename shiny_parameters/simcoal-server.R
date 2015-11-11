#
# file to be included in server.R that specifies parts of the simcoal backend
#

hst <- reactive({

    if (!is.null(input$histplotDblclick)) lstdblclick <<- input$histplotDblclick
    if (!is.null(input$histplotClick)) lstclick <<- input$histplotClick

    if (!is.null(rValues$history))
        {
            plist <- unique(c(rValues$history[,2],rValues$history[,3]))
#            if (length(plist)!=input$numpops) {
            if (length(plist)!=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops) {
                rValues$history <- NULL
            }
        }
    
    if (is.null(rValues$history))
        {
#            if (is.null(input$numpops)) {pops <- 4} else {pops <- input$numpops}
            pops <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops
            rValues$history <-create.new.history(npop=pops)
        }  else  {
            h <- rValues$history
            rValues$history <-simcoal.history.change(rValues$history,lstclick,
                                                     lstdblclick)
            if (!identical(h,rValues$history))
                {
                    lstdblclick <<- NULL
                    lstclick <<- NULL
                }
        }
    rValues$history
})

output$simhistPlot <- renderPlot({
            simcoal.history.plot(hst())
})

output$simhistTbl <- renderTable({
    df <- hst()
    rownames(df) <- 1:dim(df)[1]
    df
})


output$clickinfo <- renderText({
    c(paste("clickx",input$histplotClick$x),"\n",
      paste("dblclickx",input$histplotDblclick$x),"\n",
      paste("dblclicky",input$histplotDblclick$y))
})
