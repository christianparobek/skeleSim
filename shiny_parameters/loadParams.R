output$uiSelectParamObj <- renderUI({
  if(!is.null(input$fileParams)) {
#      print(ls(supportValues$ssLoadEnv))
      rm(list=ls(envir=supportValues$ssLoadEnv),envir=supportValues$ssLoadEnv)
      load(input$fileParams$datapath, envir = supportValues$ssLoadEnv)
#      print(ls(supportValues$ssLoadEnv))
      objs <- ls(envir = supportValues$ssLoadEnv)
      is.skeleSim <- sapply(objs, function(x) {
          class(eval(parse(text = x), envir = supportValues$ssLoadEnv)) == "skeleSim.params"
      })
      rm(list = objs[!is.skeleSim], envir = supportValues$ssLoadEnv)
      objs <- objs[is.skeleSim]
      if(length(objs) > 0) {
          obj.list <- as.list(objs)
          names(obj.list) <- objs
          selectInput(
              "slctParams",
              label = h5("Selected parameter object"),
              choices = obj.list
              )
      } else h5("<No skeleSim parameter objects found>")
  } else NULL
})

observeEvent(input$slctParams,{
  currentTitle <- if(!is.null(input$slctParams)) {
    rValues$ssClass <- get(input$slctParams, envir = supportValues$ssLoadEnv)
    rValues$ssClass@title
  } else {
    ""
  }
  updateTextInput(session, "txtTitle", value = rValues$ssClass@title)
})

