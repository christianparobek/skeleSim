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
    } else h5("<EMPTY>")
  } else h5("<EMPTY>")
})

# observe({
#   ssClass <- get(input$slctParams)
# })
