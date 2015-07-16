output$uiSelectParamObj <- renderUI({
  if(!is.null(input$fileParams)) {
    load(input$fileParams$datapath, envir = ssUserEnv)
    objs <- ls(envir = ssUserEnv)
    is.skeleSim <- sapply(objs, function(x) {
      class(eval(parse(text = x), envir = ssUserEnv)) == "skeleSim.params"
    })
    rm(list = objs[!is.skeleSim], envir = ssUserEnv)
    objs <- objs[is.skeleSim]
    if(length(objs) > 0) {
      obj.list <- as.list(objs)
      names(obj.list) <- objs
      selectInput(
        "slctParams",
        label = h5("Select parameter object"),
        choices = obj.list
      )
    } else h5("<No skeleSim parameter objects found>")
  } else NULL
})

observe({
  if(!is.null(input$slctParams)) {
    ssClass <<- get(input$slctParams, envir = ssUserEnv)
    output$txtChangeTitle <- renderUI({
      textInput("txtNewTitle", label = h4("Title"), value = ssClass@title)
    })
    print(ssClass@title)
  }
})

