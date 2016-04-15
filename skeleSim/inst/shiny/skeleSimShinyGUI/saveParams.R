###shinyFiles approach to saving the current skelesim object
shinyFileSave(input, 'ssClassSave', filetypes=c("rdata","RData","RDATA","rData","rda"), session=session, roots=VolumeRoots)
###

observeEvent(input$ssClassSave,
             {
                 if (debug()) print("parsing")
                 path <- as.character(parseSavePath(VolumeRoots,input$ssClassSave)$datapath)
                 if (debug()) print(path)
                 if (debug()) print(normalizePath(path))
                 if (debug()) print("done parsing")
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




        
observeEvent(input$btnSave,{
    req(rValues$ssClass)
    req(supportValues$simroot)
    d <- getwd()
    setwd(supportValues$simroot)
    if (rValues$ssClass@simulator=="fsc")
        {
            if (debug()) print("running fsc.write")
            fsc.write(rValues$ssClass)
        }

    if (rValues$ssClass@simulator=="rms")
    {
        rms.write(rValues$ssClass)
    }
    setwd(d)
})
