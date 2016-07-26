###
### implements a lot of the backend for the viz tab
###

###  selections set up when object read in

output$globalMainPlot <- renderPlot({
    req(rValues$ssClass@analysis.results)
    if ("Global" %in% names(rValues$ssClass@analysis.results[[1]]))
        if (input$scompare=="s")
            gg.global(rValues$ssClass,stats=input$gstatsel,scenario=input$vizScenario)
        else
            gg.global.scmp(rValues$ssClass,stats=input$gstatsel)
})

output$locusMainPlot <- renderPlot({
    req(rValues$ssClass@analysis.results)
    if ("Locus" %in% names(rValues$ssClass@analysis.results[[1]]))
          if (input$scompare=="s")
            gg.locus(rValues$ssClass,stats=input$lstatsel,scenario=input$vizScenario)
        else
            gg.locus.scmp(rValues$ssClass,stats=input$lstatsel)
})

output$pairwiseMainPlot <- renderPlot({
    req(rValues$ssClass@analysis.results)
    if ("Pairwise" %in% names(rValues$ssClass@analysis.results[[1]]))
          if (input$scompare=="s")
          {
              gg.pairwise(rValues$ssClass,stats=input$pstatsel,scenario=input$vizScenario)
          } else {
              gg.pairwise.scmp(rValues$ssClass,stats=input$pstatsel)
          }
})


output$dlGlobal <- downloadHandler(
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
        paste0(rValues$ssClass@title,"-Global-",gsub(" |:","_",date()),".csv")
    },
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      write.table(df.global(rValues$ssClass), file, sep = ",",
        row.names = FALSE)
    },
    contentType="text/csv"
  )

output$dlLocus <- downloadHandler(
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
        paste0(rValues$ssClass@title,"-Locus-",gsub(" |:","_",date()),".csv")
    },
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      write.table(df.locus(rValues$ssClass), file, sep = ",",
        row.names = FALSE)
    },
    contentType="text/csv"

  )

output$dlPairwise <- downloadHandler(
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
        paste0(rValues$ssClass@title,"-Pairwise-",gsub(" |:","_",date()),".csv")
    },
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      write.table(df.pairwise(rValues$ssClass), file, sep = ",",
        row.names = FALSE)
    },
    contentType="text/csv"

  )


