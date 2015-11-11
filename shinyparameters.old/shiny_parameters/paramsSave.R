output$uiBtnSaveParams <- renderUI({
  # sets global variable 'objLabel' based on user update of class title
  objLabel <<- if(is.null(input$txtTitle) | input$txtTitle == "") {
    NULL
  } else {
    ssClass@title <<- input$txtTitle
    make.names(input$txtTitle)
  }

  # reset sim running button
  output$btnRun <- renderUI({h5("<Parameter File Not Saved>")})

  # fill construct save button and return label field
  if(is.null(objLabel)) {
    h5("<EMPTY>")
  } else {
    actionButton("btnSaveParams", objLabel)
  }
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

    # report write status and setup sim running button
    output$txtSaveStatus <- renderText({
      if(file.exists(paramsFname)) {
        output$btnRun <- renderUI({
          actionButton("btnRunSim", paste("Run", fnameLabel))
        })
        paste("Wrote skeleSim parameter file: '", paramsFname, "'", sep = "")
      } else "Error writing files."
    })
  }
})

