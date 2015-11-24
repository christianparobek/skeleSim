output$txtObjLabel <- renderText({
  objLabel <<- if(is.null(input$txtTitle) | input$txtTitle == "") {
    NULL
  } else {
    make.names(input$txtTitle)
  }
  if(is.null(objLabel)) "<EMPTY>" else objLabel
})

observeEvent(input$btnSaveParams, {
  # uses global variables 'objLabel' and 'fnameLabel' defined in 'global.R'
  if(!is.null(objLabel)) {
    # create 'fnameLabel' and parameter filename (paramsFname)
    fnameLabel <<- paste(objLabel, format(Sys.time(), "%Y%m%d.%H%M"), sep = ".")
    paramsFname <- paste(fnameLabel, "params.rdata", sep = ".")
    # create parameter object based on user-defined 'objLabel' and save to file
    assign(objLabel, ssClass)
    save(list = objLabel, file = paramsFname)

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
  }
})