###shinyFiles approach to saving the current skelesim object
shinyFileSave(input, 'ssClassSave', filetypes=c("rdata","RData","RDATA","rData","rda"), session=session, roots=VolumeRoots)
###

observeEvent(input$ssClassSave,
             {
                 if (debug()) print("parsing")
                 path <- as.character(parseSavePath(VolumeRoots,input$ssClassSave)$datapath)
                 fn <- basename(path)
                 dir <- normalizePath(dirname(path))
                 if (debug()) print(path)
                 if (debug()) print(normalizePath(path))
                 if (debug()) print("done parsing")
                 assign("ssClass",rValues$ssClass)
                 save(file=paste0(dir,"/",fn),list="ssClass")
             })

output$txtObjLabel <- renderText({
  supportValues$objLabel <- if(is.null(req(rValues$ssClass@title)) | rValues$ssClass@title == "") {
    NULL
  } else {
    make.names(rValues$ssClass@title)
  }
  if(is.null(supportValues$objLabel)) ret <- "<NONE>" else ret <- supportValues$objLabel
  paste("Title of skelesim object currently in memory:",ret)
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

observeEvent(input$resetButton,{
#    rValues <- reactiveValues(ssClass=ssClassInit(),
#                          scenarioNumber=1,
#                          lstScenario=1,
#                          migrationNumber=0,
#                          lstMigration = 0,
#                          localDemoNumber = 1,
#                          lstLocalDemo = 1,
#                          EpochNumber = 1,
#                          lstEpoch = 1,
#                          history=NULL,
#                          msg=NULL)

    rValues$ssClass <- ssClassInit()
    rValues$scenarioNumber <- 1
    rValues$lstScenario <- 1
    rValues$migrationNumber=0
    rValues$lstMigration = 0
    rValues$localDemoNumber = 1
    rValues$lstLocalDemo = 1
    rValues$EpochNumber = 1
    rValues$lstEpoch = 1
    rValues$history=NULL
    rValues$msg=NULL
                          
    updateUIs()
})
