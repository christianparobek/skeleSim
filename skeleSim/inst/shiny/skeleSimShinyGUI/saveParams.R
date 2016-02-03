###shinyFiles approach to saving the current skelesim object
shinyFileSave(input, 'ssClassSave', filetypes=c("rdata","RData","RDATA","rData","rda"), session=session, roots=VolumeRoots)
###

observeEvent(input$ssClassSave,
             {
                 print("parsing")
                 path <- as.character(parseSavePath(VolumeRoots,input$ssClassSave)$datapath)
                 print(path)
                 print(normalizePath(path))
                 print("done parsing")
                 assign("ssClass",rValues$ssClass)
                 save(file=normalizePath(path),list="ssClass")
             })

output$txtObjLabel <- renderText({
  supportValues$objLabel <- if(is.null(req(rValues$ssClass@title)) | rValues$ssClass@title == "") {
    NULL
  } else {
    make.names(rValues$ssClass@title)
  }
  if(is.null(supportValues$objLabel)) ret <- "<NONE>" else ret <- supportValues$objLabel
  paste("Title of skelesim object running:",ret)
})


