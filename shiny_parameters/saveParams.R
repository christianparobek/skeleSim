

output$txtObjLabel <- renderText({
  supportValues$objLabel <- if(is.null(rValues$ssClass@title) | rValues$ssClass@title == "") {
    NULL
  } else {
    make.names(rValues$ssClass@title)
  }
  if(is.null(supportValues$objLabel)) ret <- "<NONE>" else ret <- supportValues$objLabel
  paste("Current title of skelesim object:",ret)
})

observeEvent(input$btnSaveParams, {

    if(!is.null(supportValues$objLabel)) {
        ## create 'fnameLabel' and parameter filename (paramsFname)
        supportValues$fnameLabel <- paste(supportValues$objLabel, format(Sys.time(), "%Y%m%d.%H%M"), sep = ".")
        paramsFname <- paste(supportValues$fnameLabel, "params.rdata", sep = ".")
        paramsFname <- paste0("../",paramsFname)
        print(paramsFname)
        ## create parameter object based on user-defined 'objLabel' and save to file

        assign("ssClass", rValues$ssClass)
        save(list = "ssClass" , file = paramsFname) #this needs updating to be resiliant

        
        output$txtSaveStatus <- renderText({
            if(file.exists(paramsFname)) {
                paste("Wrote skeleSim parameter file: '", paramsFname, "'", sep = "")
            } else {
                "Error writing files."
            }
        })

        output$btnRun <- renderUI({
            actionButton("btnRunSim", paste("Run", supportValues$fnameLabel))
        })
    }
})
