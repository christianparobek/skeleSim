  megalistValues <- reactive({
    # Compose data frame
    data.frame(
      Name = c(#"SNPs", 
               "Non-Diploid",
               "Many Markers",
               "Large Pop Size",
               "Complex History",
               "Deep Timeframes",
               "Model Demography",
               "Management Question",
               "Fast Completion",
               "Large Computer"),
      Value = as.character(c(#input$snps, 
                             input$non.diploid,
                             input$marker.num,
                             input$pop.size,
                             input$complex.hist,
                             input$deep.time,
                             input$demography,
                             input$management,
                             input$completion.time,
                             input$computer)), 
      stringsAsFactors=FALSE)
  })

  # Server-side support for rendering megalist values
  output$values <- renderTable({
    megalistValues()
  })
  
  # Server-side support for rendering simulation name
  output$simname <- renderText({
    input$simname
  })
  
  # Server-side support for calculating coal / foward score
  
#  output$simscore <- 
    
    # default responses
#    responses <- c(input$snps, 
#                   input$non.diploid,
#                   input$marker.num,
#                   input$pop.size,
#                   input$complex.hist,
#                   input$deep.time,
#                   input$demography,
#                   input$management,
#                   input$completion.time,
#                   input$computer)


  # response weights
#  forward.wts <- c(0, 0, 0.3, 0.2, 0, 0.2, 1, 1, 0.2, 0.3)
  
  # get relative 'score' for each model
#  fwd.score <- sum(forward.wts * responses) / length(responses)
#  cat("Coalescent score: ", 1 - fwd.score, "\n")
#  cat("Forward-time score: ", fwd.score, "\n")
