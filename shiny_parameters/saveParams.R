observeEvent(input$btnSaveParams, {
  # save skeleSim class object (ssClass)
  #   uses global variable 'fnameLabel' defined in 'global.R'
  fnameLabel <<- gsub("[[:punct:]|[[:space:]]", ".", input$txtParamsFname)
  fnameLabel <<- paste(fnameLabel, format(Sys.time(), "%Y%m%d.%H%M"), sep = ".")
  paramsFname <- paste(fnameLabel, ".params.rdata", sep = "")
# !!! CHANGE 'label' to ssClass !!!
  save(fnameLabel, file = paramsFname)

  output$txtSaveStatus <- renderText({
    if(file.exists(paramsFname)) {
      paste("Wrote skeleSim parameter file: '", paramsFname, "'", sep = "")
    } else {
      "Error writing files."
    }
  })

  output$btnRun <- renderUI({
    actionButton("btnRunSim", paste("Run", fnameLabel))
  })
})