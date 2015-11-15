observeEvent(input$btnRunSim, {
  # write script
  #   uses global variable 'fnameLabel' defined in 'global.R'
  if(!is.null(supportValues$fnameLabel)) {
    scriptFname <- paste(supportValues$fnameLabel, ".script.R", sep = "")
    paramFname <- paste(supportValues$fnameLabel, ".params.rdata", sep = "")
    outFname <- paste(supportValues$fnameLabel, ".skeleSim.out.rdata", sep = "")
    
    write("setwd('../')", file = scriptFname)
    write("getwd()",file = scriptFname, append=TRUE)
    write("rm(list = ls())", file = scriptFname, append=TRUE)
    write("library(methods)", file = scriptFname, append = TRUE)

    srcfiles <- ifelse(rValues$ssClass@simulator=="fsc",'source("fastsimcoal.skeleSim.R")','source("rmetasim.skeleSim.R")')
    write('source("skeleSim.classes.R")', file = scriptFname, append = TRUE)
    write('source("skeleSim.funcs.R")', file = scriptFname, append = TRUE)
    write(srcfiles, file = scriptFname, append = TRUE)
print("srcfiles")
    print(srcfiles)
    print("printed srcfiles")
    line <- paste("load('", paramFname, "')", sep = "")
    write(line, file = scriptFname, append = TRUE)

    write("ls()", file = scriptFname, append = TRUE)
    
    write("ssClass <- runSim(ssClass)", file = scriptFname, append = TRUE)
#    write("ssClass@timing <- NULL", file = scriptFname, append = TRUE)
    line <- paste("save(ssClass, file = '", outFname, "')", sep = "")
    write(line, file = scriptFname, append = TRUE)

    cond <- TRUE
    output$txtRunStatus <- renderText({
      if(cond) {
        cmd <- paste("nohup R CMD BATCH --no-restore", scriptFname)
        system(cmd, wait = FALSE)
        "Simulation run commenced."
        supportValues$fnameLabel <- NULL
        stopApp()
      } else {
        "Error starting simulations."
      }
    })

  }
})
