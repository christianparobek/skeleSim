
###choose working directory using the shinyFiles package
shinyDirChoose(input,"workFolder",session=session,roots=VolumeRoots,hidden=F)
###
####choose file to load using the shinyFiles package
#roots <- getVolumes()
#observe
shinyFileChoose(input,"fileParams",filetypes=c("rdata","RData","RDATA","rData","rda"),session=session,roots=VolumeRoots)
#done with the shinyFiles file selection code

###make sure that changing the working directory updates the supportValues$simroot slot, but not any other time
observeEvent(input$workFolder,{
  #  if (!is.null(input$workFolder))
        supportValues$simroot <- normalizePath(gsub("/+","/",parseDirPath(VolumeRoots,input$workFolder)))
})
               
output$uiSelectParamObj <- renderUI({
  if(!is.null(input$fileParams)) {
      fileParams <- parseFilePaths(VolumeRoots,input$fileParams) #convert shinyFiles object into more familiar inputFiles format
#      print(fileParams$datapath)
      rm(list=ls(envir=supportValues$ssLoadEnv),envir=supportValues$ssLoadEnv)
      load(file=as.character(fileParams$datapath), envir = supportValues$ssLoadEnv)
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
              label = h5("Selected parameter object (could be multiple objects per file)"),
              choices = isolate(obj.list)
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
  updateTextInput(session, "txtTitle", value = isolate(rValues$ssClass@title))
  updateUIs()
})

