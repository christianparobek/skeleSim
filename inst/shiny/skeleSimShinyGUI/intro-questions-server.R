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
        Value = c(#input$snps, 
            input$non.diploid,
            input$marker.num,
            input$pop.size,
            input$complex.hist,
            input$deep.time,
            input$demography,
            input$management,
            input$completion.time,
            input$computer), 
        stringsAsFactors=FALSE)
  })

  # Server-side support for rendering megalist values
  output$values <- renderTable({
    megalistValues()
  })
  
 # # Server-side support for rendering simulation name
 # output$simname <- renderText({
 #   input$simname
 # })
  
  # Server-side support for calculating coal / foward score

simscore <- reactive(
    {
        df <- megalistValues()
#        forward.wts <- c(0.5, 0.5, 0.2, 0.2, 0.2, 1, 1, 0.2, 1)

#        forward.wts <- rep(1,9)
        forward.wts <- c(0.5, 0.5, 0.2, 0.5, 0.1, 1, 1, 0.1, 1 )

        fwd <- df[,2]*forward.wts

        round(mean(fwd),3)
    })

output$coalscore <- renderText({paste("Coalescent weighting score",1-simscore())})
output$forwardscore <- renderText({paste("Forward-time weighting score",simscore())})

observeEvent(simscore(),{
    updateCheckboxInput(session,"coalescent","Coalescent simulator?",value=simscore()<0.5)
})
