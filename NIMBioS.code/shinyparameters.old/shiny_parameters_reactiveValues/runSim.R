observeEvent(input$btnRunSim, {
  # write script
  #   uses global variable 'fnameLabel' defined in 'global.R'
  if(!is.null(fnameLabel)) {
    scriptFname <- paste(fnameLabel, ".script.R", sep = "")
    paramFname <- paste(fnameLabel, ".rdata", sep = "")
    outFname <- paste(fnameLabel, ".skeleSim.out.rdata", sep = "")
    write("rm(list = ls())", file = scriptFname)
    write("library(methods)", file = scriptFname, append = TRUE)
    paramFname <- "test.ssClass.rdata"
    line <- paste("load('", paramFname, "')", sep = "")
    write(line, file = scriptFname, append = TRUE)
    write("ssClass <- runSim(ssClass)", file = scriptFname, append = TRUE)
    write("ssClass@timing <- NULL", file = scriptFname, append = TRUE)
    line <- paste("save(ssClass, file = '", outFname, "')", sep = "")
    write(line, file = scriptFname, append = TRUE)

    cond <- TRUE
    output$txtRunStatus <- renderText({
      if(cond) {
        cmd <- paste("nohup R CMD BATCH --no-restore", scriptFname)
        system(cmd, wait = FALSE)
        "Simulation run commenced."
        fnameLabel <<- NULL
        stopApp()
      } else {
        "Error starting simulations."
      }
    })

  }
})