#
# file to be included in server.R that specifies parts of the simcoal backend
#

hst <- reactive({

    if (!is.null(input$histplotDblclick)) lstdblclick <<- input$histplotDblclick
    if (!is.null(input$histplotClick)) lstclick <<- input$histplotClick

    if (!is.null(rValues$history))
        {
            plist <- unique(c(rValues$history[,2],rValues$history[,3]))
#            if (length(plist)!=input$numpops) {
            if (length(plist)!=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops) {
                rValues$history <- NULL
            }
        }
    
    if (is.null(rValues$history))
        {
#            if (is.null(input$numpops)) {pops <- 4} else {pops <- input$numpops}
            pops <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops
            rValues$history <-create.new.history(npop=pops)
        }  else  {
            h <- rValues$history
            rValues$history <-simcoal.history.change(rValues$history,lstclick,
                                                     lstdblclick)
            if (!identical(h,rValues$history))
                {
                    lstdblclick <<- NULL
                    lstclick <<- NULL
                }
        }
    rValues$history
})

output$simhistPlot <- renderPlot({
            simcoal.history.plot(hst())
})

output$simhistTbl <- renderTable({
    df <- hst()
    rownames(df) <- 1:dim(df)[1]
    df
})


output$clickinfo <- renderText({
    c(paste("clickx",input$histplotClick$x),"\n",
      paste("dblclickx",input$histplotDblclick$x),"\n",
      paste("dblclicky",input$histplotDblclick$y))
})


samp.times <- function()
    {
        print("running samptime")
        if ((is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@sample.times))|
            (length(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@sample.times)!=
             rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops))
            ret <- rep(0,rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops)
        else 
            {
                ret <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@sample.times
            }
        matrix(ret,nrow=1)
    }

output$samptime <- renderUI({
    print("creating st vector")
    matrixInput("stvec","Vector of sampling times (corresponds to populations) (please don't use +/- buttons at right [temporary])",
                as.data.frame(samp.times()))
})

growth.rates <- function()
    {
        print("running growth.rates")
        if ((is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@growth.rate))|
            (length(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@growth.rate)!=
             rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops))
            ret <- rep(0,rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops)
        else 
            {
                ret <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@growth.rate
            }
        matrix(ret,nrow=1)
    }

output$growthrate <- renderUI({
    print("creating growthrate vector")
    matrixInput("grvec","Vector of growth rates (corresponds to populations) (please don't use +/- buttons at right [temporary])",
                as.data.frame(growth.rates()))
})


output$simexec <- renderUI({
    ui <- textInput("fscexec", "No fastsimcoal executable in path: enter value", value = "")
    if (!is.null(rValues$ssClass@simulator))
        if (rValues$ssClass@simulator=="fsc")
            if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params))
                    {
                        sim.exec <- c(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@fastsimcoal.exec,
                                      supportValues$simexec)
                        sim.exec <- sim.exec[!is.null(sim.exec)]
                        sim.exec <- unique(sim.exec)
                        sim.exec <- basename(Sys.which(sim.exec))
                        sim.exec <- sim.exec[nchar(sim.exec)>0]
                        ui <- selectInput("fscexec","Select simcoal executable",choices=sim.exec)
                    }
    ui
})


###########locus params (actually derived from scenario-specific information)
observeEvent(rValues$ssClass@scenarios[[rValues$scenarioNumber]],{
    if (!is.null(rValues$ssClass@simulator))
        if (rValues$ssClass@simulator=="fsc")
            {
                nl <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci
                print (paste("numloci",nl))
                mat <- matrix("",nrow=nl,ncol=5)
                for (l in 1:nl)
                    {
                        mat[l,1] <- ifelse(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@locus.type=="sequence","DNA","MICROSAT")
                        mat[l,2] <- ifelse(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@locus.type=="sequence",
                                           as.character(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@sequence.length),1)
                        mat[l,3] <- as.character(0)
                        mat[l,4] <- as.character(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mut.rate[l])
                        mat[l,5] <- as.character(1/3)
                    }
                rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@locus.params <- mat
            }
})
