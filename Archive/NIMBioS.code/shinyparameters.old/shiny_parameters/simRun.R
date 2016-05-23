observeEvent(input$btnRunSim, {
  # write script
  #   uses global variable 'fnameLabel' defined in 'global.R'
  if(!is.null(fnameLabel)) {
    scriptFname <- paste(fnameLabel, ".script.R", sep = "")
    paramFname <- paste(fnameLabel, ".rdata", sep = "")
    resultFname <- paste(fnameLabel, ".result.rdata", sep = "")
    write("rm(list = ls())", file = scriptFname)
    #write("library(skeleSim)", file = scriptFname, append = TRUE)
    write("library(methods)", file = scriptFname, append = TRUE)
    line <- paste("load('", paramFname, "')", sep = "")
    write(line, file = scriptFname, append = TRUE)
    result <- paste(objLabel, ".result", sep = "")
    line <- paste(result, " <- runSim(", objLabel, ")", sep = "")
    write(line, file = scriptFname, append = TRUE)
    line <- paste("save(", result, ", file = '", resultFname, "')", sep = "")
    write(line, file = scriptFname, append = TRUE)

    cond <- TRUE
    output$txtRunStatus <- renderText({
      if(cond) {
        cmd <- paste("nohup R CMD BATCH --no-restore", resultFname)
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