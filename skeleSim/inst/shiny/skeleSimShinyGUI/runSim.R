observeEvent(input$btnRun, {
    print("in run")
    if(!is.null(supportValues$objLabel)) {
        print("past first null test")
        print("supportValues$objLabel")
        ## create 'fnameLabel' and parameter filename (paramsFname)
        supportValues$fnameLabel <- paste(supportValues$objLabel, format(Sys.time(), "%Y%m%d.%H%M"),round(runif(1,min=0,max=10000)), sep = ".")
        paramsFname <- paste(supportValues$fnameLabel, "params.rdata", sep = ".")
        paramsFname <- paste0(ifelse(is.null(supportValues$simroot),".",supportValues$simroot),"/",paramsFname)
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

#        output$btnRun <- renderUI({
#            actionButton("btnRunSim", paste("Run", supportValues$fnameLabel))
#        })

#### do the running now
        
        if(!is.null(supportValues$fnameLabel)) {
            scriptFname <- paste(supportValues$fnameLabel, ".script.R", sep = "")
            paramFname <- paste(supportValues$fnameLabel, ".params.rdata", sep = "")
            outFname <- paste(supportValues$fnameLabel, ".skeleSim.out.rdata", sep = "")
            
            scriptFname <- paste(supportValues$simroot,scriptFname,sep="/")
            
            write("rm(list = ls())", file = scriptFname, append=TRUE)
            write("library(methods)", file = scriptFname, append = TRUE)
            
####these next lines will need to be modified once this code is made into a package
### right now we are depending on the fact that shiny apps run in the directory where
### server.R is located

### Right now the way I'm going about things is to switch to the skeleSim root directory
#### source all the relevant files, and then change to the directory where the simulation
####  should be run (supportValues$simroot)

            write(paste0("library(skeleSim)"), file = scriptFname)
            write("getwd()",file = scriptFname, append=TRUE)
            
#            srcfiles <- ifelse(rValues$ssClass@simulator=="fsc",
#                               'source("fastsimcoal.skeleSim.R")',
#                               'source("rmetasim.skeleSim.R")'
#                               )
#            write('source("skeleSim.classes.R")', file = scriptFname, append = TRUE)
#            write('source("skeleSim.funcs.R")', file = scriptFname, append = TRUE)
#            write(srcfiles, file = scriptFname, append = TRUE)
#            print("srcfiles")
#            print(srcfiles)
#            print("printed srcfiles")
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
            print("about to run new R session")
            print(ifelse(supportValues$OS=="unix","Unix-based system","Windows-based system"))
            if (supportValues$OS=="unix") #
                {
                    cmd <- paste("nohup R CMD BATCH --no-restore", scriptFname)
                    system(cmd, wait = FALSE)
                    
                } else {
                    cmd <- paste("start R CMD BATCH --no-restore", scriptFname)
                    shell(cmd,wait=F)
                }
            print("about to run system")

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


