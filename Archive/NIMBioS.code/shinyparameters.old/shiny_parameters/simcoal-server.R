#
# file to be included in server.R that specifies parts of the simcoal backend
#

hst <- reactive({
    if (!exists("histry"))
        {
            histry <<- NULL
        } 

    if (!is.null(input$histplotDblclick)) lstdblclick <<- input$histplotDblclick
    if (!is.null(input$histplotClick)) lstclick <<- input$histplotClick

    if (!is.null(histry))
        {
            plist <- unique(c(histry[,2],histry[,3]))
            if (length(plist)!=input$numpops) {
#                print(plist)
#                print(input$numpops)
#                print ("Resetting history")
                histry <<- NULL
            }
        }
    

    if (is.null(histry))
        {
            if (is.null(input$numpops)) {pops <- 4} else {pops <- input$numpops}
            histry<<-create.new.history(npop=pops)
#            print("TEST")
        }  else  {
            h <- histry
            histry<<-simcoal.history.change(histry,lstclick,
                                            lstdblclick)
            if (!identical(h,histry))
                {
                    lstdblclick <<- NULL
                    lstclick <<- NULL
                }
        }
    histry
})

output$simhistPlot <- renderPlot({
            simcoal.history.plot(hst())
})

output$clickinfo <- renderText({
    c(paste("clickx",input$histplotClick$x),"\n",
      paste("dblclickx",input$histplotDblclick$x),"\n",
      paste("dblclicky",input$histplotDblclick$y))
})

output$simhistTbl <- renderTable({
    df <- hst()
    rownames(df) <- 1:dim(df)[1]
    df
})
