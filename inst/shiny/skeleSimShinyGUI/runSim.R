output$runButton <- renderUI({
    if (!is.null(req(supportValues$simroot)))
    actionButton("btnRun","Run simulation")
})
output$saveButton <- renderUI({
    if (!is.null(req(supportValues$simroot)))
    actionButton("btnSave","Save example inputs for each scenario")
})

output$rootText <- renderText({paste("Simulation inputs, outputs and work directory will be written in:", supportValues$simroot)})



observeEvent(input$btnRun, {
    req(supportValues$simroot)
#    print("in run")
    if(!is.null(supportValues$objLabel)) {
#        print("past first null test")
#        print("supportValues$objLabel")
        ## create 'fnameLabel' and parameter filename (paramsFname)
        supportValues$fnameLabel <- paste(supportValues$objLabel, format(Sys.time(), "%Y%m%d.%H%M"),round(runif(1,min=0,max=10000)), sep = ".")
        paramsFname <- paste(supportValues$fnameLabel, "params.rdata", sep = ".")
        paramsFname <- paste0(ifelse(is.null(supportValues$simroot),".",supportValues$simroot),"/",paramsFname)
        if (debug()) print(paramsFname)
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

#### do the running now
        
        if(!is.null(supportValues$fnameLabel)) {
            scriptFname <- paste(supportValues$fnameLabel, ".script.R", sep = "")
            paramFname <- paste(supportValues$fnameLabel, ".params.rdata", sep = "")
            outFname <- paste(supportValues$fnameLabel, ".skeleSim.out.rdata", sep = "")
            logFname <- paste(supportValues$fnameLabel, ".log", sep = "")
            scriptFname <- paste(supportValues$simroot,scriptFname,sep="/")
            
            write("rm(list = ls())", file = scriptFname, append=TRUE)
            write("library(methods)", file = scriptFname, append = TRUE)
            

### Right now the way I'm going about things is to switch to the skeleSim root directory
#### source all the relevant files, and then change to the directory where the simulation
####  should be run (supportValues$simroot)

            write(paste0("library(skeleSim)"), file = scriptFname)
            write(paste0("library(adegenet)"),file = scriptFname, append=TRUE)
            write("getwd()",file = scriptFname, append=TRUE)
            
            cdcmd <- chartr("\\","/",paste0("setwd('",supportValues$simroot,"')")) #make this work on windows also
            write(cdcmd, file = scriptFname, append=TRUE)
            write("getwd()",file = scriptFname, append=TRUE)

#### end of the sourceing
            
            line <- paste("load('", paramFname, "')", sep = "")
            write(line, file = scriptFname, append = TRUE)
            write("ls()", file = scriptFname, append = TRUE)
            write("ssClass <- runSim(ssClass)", file = scriptFname, append = TRUE)
                                        #    write("ssClass@timing <- NULL", file = scriptFname, append = TRUE)
            line <- paste("save(ssClass, file = '", outFname, "')", sep = "")
            write(line, file = scriptFname, append = TRUE)
            
            cond <- TRUE
            if (debug()) print("about to run new R session")
            if (debug()) print(ifelse(supportValues$OS=="unix","Unix-based system","Windows-based system"))
            if (supportValues$OS=="unix") #
                {
                    cmd <- paste("cd",supportValues$simroot,"; nohup R CMD BATCH --no-restore ", scriptFname)
                    system(cmd, wait = FALSE)
                    
                } else { #windows output put in shiny directory
                    cmd <- paste("start R CMD BATCH --no-restore", scriptFname)
                    shell(cmd,wait=F)
                }
            if (debug()) print("about to run system")

####            output$txtRunStatus <- renderText({
#                if(cond) {
#                    cmd <- paste("nohup R CMD BATCH --no-restore", scriptFname)
#                    print("about to run system")
#                    system(cmd, wait = FALSE)
#                    "Simulation run commenced."
#                    supportValues$fnameLabel <- NULL
#                                        #                     stopApp()
#                } else {
#                    "Error starting simulations."
#                }
#            })
            
          }
        
        
    }
})


