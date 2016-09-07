#
# file to be included in server.R that specifies parts of the simcoal backend
#

hst <- reactive({
    if (rValues$ssClass@simulator=="fsc")
    {
        history <- NULL
        
        if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@hist.ev))
        {
            history <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@hist.ev
 #           print(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@hist.ev)
        }
        
        if (!is.null(history))
        {
            plist <- unique(c(history[,2],history[,3]))
            if (length(plist)!=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops) {
                history <- NULL
            }
 #           print("past plist")
if (debug())              print(history)
        }
        
        if (is.null(history))
        {
            pops <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops
            history <-create.new.history(npop=pops)
        }  else  {
            h <- history
#            print("about to change")
            if (!is.null(pointValues$dblclick))
                if (!is.null(pointValues$click))
                {
#                    print("altering history")
if (debug())                      print(pointValues$click$x)
if (debug())                      print(paste(pointValues$dblclick$x,pointValues$dblclick$y))
                    history <-simcoal.history.change(history,pointValues$click,
                                                     pointValues$dblclick)
if (debug())                      print(history)
                }
#            print("just ran change")
            
            pointValues$click <- NULL
            pointValues$dblclick <- NULL
            
        }
#        print("returning from hst()")
        history
    }
})

output$simhistPlot <- renderPlot({
    if (debug()) print("about to plot history")
    h <- hst()
    if (debug()) print(h)
    simcoal.history.plot(h)
})

output$simhistTbl <- renderTable({
    df <- hst()
    rownames(df) <- 1:dim(df)[1]
    df
})

output$simhistEditTbl <- renderUI({
#    print("creating simhistEditTbl")
    matrixInput("simhist","time | source | sink | migrants | new.size | growth.rate | migr.matrix",
                as.data.frame(hst()))
})

observeEvent(input$histplotClick,
{
    pointValues$click <- input$histplotClick
})

observeEvent(input$histplotDblClick,
{
    pointValues$dblclick <- input$histplotDblClick
    if (!is.null(pointValues$click))
        if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@hist.ev))
            {
                if (debug())  print("inside observe event dblclick, inside all not nulls")
                h=hst()
                if (!historiesEqual(h,
                                    rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@hist.ev))
                {
                    if (debug())  print("histories are not equal")
                    rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@hist.ev <- h
                }
            }

})

                                        #
#observeEvent(pointValues,{
#    print("in observevent pointValues")
#    if (!is.null(pointValues$click))
#            if (!is.null(pointValues$dblclick))
#                if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@hist.ev))
#            {
#                h=hst()
#                if (!historiesEqual(h,
#                                   rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@hist.ev))
#                {
#                    rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@hist.ev <- h
#                }
#            }
#    
#})


observeEvent(input$simhist,{
    if (debug()) print("input$simhist modified")
    simhist <- input$simhist
    mnum <- 0
    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration))
        mnum <- length(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration)
    if (debug()) print("assigned mnum")
    ps <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@pop.size
    if (debug()) print(paste("got popsize",paste(ps,collapse=",")))
    if (!isTRUE(all.equal(req(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@hist.ev),
                          simhist)))
    {
        if (debug()) print("hist modified to new value")
#        print("this is the hist.ev")
        if (dim(simhist)[1]==0) simhist <- NULL
#        print(simhist)
#        print(dim(simhist))
        hevck <- fsc.histEvCheck(hist.ev=simhist,
                            pop.size=ps,
#                            0,
                            growth.rate=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@growth.rate,
                            num.mig.mats=length(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration))
if (debug())          if (length(hevck)==0) print ("hevck not set") else print(paste("hevck",hevck))
        if ((hevck))
            rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@hist.ev  <- simhist
        else
            output$simhistEditTbl <- renderUI({ #redraw matrix with stored values, the input$simhist values are not legal
                matrixInput("simhist","time | source | sink | migrants | new.size | growth.rate | migr.matrix",
                            as.data.frame(hst()))
            })

    }
})

observeEvent(input$addHistEvent,{
    rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@hist.ev <-
        rbind(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@hist.ev,
              c(1,0,0,1,0,0,0))
})

observeEvent(input$removeLastEvent,{
    if (debug())  print("in observEvent removeLastEvent")
    hist <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@hist.ev
    if (fsc.histEvCheck(hist[-1,],
                   rValues$ssClass@scenarios[[rValues$scenarioNumber]]@pop.size,
                   rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@growth.rate,
                   num.mig.mats=length(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration)))
        if ((dim(hist)[2]-1)>=max(c(hist[,2:3]))) #if the dimensions of the matrix are large enough for every pop to coalesce
        {
            if (debug())  print("changing hsit.ev as a consequence of removing a row")
            if (debug())  print(hist)
            rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@hist.ev <- hist[-dim(hist)[1],]
        }
})


output$clickinfo <- renderText({
    c(paste("clickx",input$histplotClick$x),"\n",
      paste("dblclickx",input$histplotDblclick$x),"\n",
      paste("dblclicky",input$histplotDblclick$y))
})


samp.times <- function()
    {
        if (debug()) print("running samptime")
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
#    print("creating st vector")
    matrixInput("stvec","Vector of sampling times (corresponds to populations)",
                as.data.frame(samp.times()))
})

growth.rates <- function()
    {
        if (debug()) print("running growth.rates")
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
#    print("creating growthrate vector")
    matrixInput("grvec","Vector of growth rates (corresponds to populations)",
                as.data.frame(growth.rates()))
})


#output$simexec <- renderUI({
#    ui <- textInput("fscexec", "No fastsimcoal executable in path: enter value", value = "")
#    if (!is.null(rValues$ssClass@simulator))
#        if (rValues$ssClass@simulator=="fsc")
#            if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params))
#                    {
#                        sim.exec <-
#                            c(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@fastsimcoal.exec,
#                                      supportValues$simexec)
#                        sim.exec <- sim.exec[!is.null(sim.exec)]
#                        sim.exec <- basename(Sys.which(sim.exec))
#                        sim.exec <- unique(sim.exec)
#                        sim.exec <- sim.exec[nchar(sim.exec)>2]
#                        if (debug()) print("rendering UI for simexec")
#                        if (debug()) print(sim.exec)
#                        ui <- selectInput("fscexec","Select fastsimcoal executable",selected=sim.exec[1],
#                                          choices=sim.exec)
#                    }
#    ui
#})


###########locus params (actually derived from scenario-specific information)
observeEvent(rValues$ssClass@scenarios[[rValues$scenarioNumber]],{
    if (!is.null(rValues$ssClass@simulator))
        if (rValues$ssClass@simulator=="fsc")
            {
                nl <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci
                if (debug()) print (paste("numloci",nl))
                mat <- matrix("",nrow=nl,ncol=5)
                for (l in 1:nl)
                    {
                        mat[l,1] <- ifelse(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@locus.type%in%c("sequence","SNP"),"DNA","MICROSAT")
                        mat[l,2] <- ifelse(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@locus.type=="sequence",
                                           as.character(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@sequence.length),1)
                        mat[l,3] <- as.character(0)
                        mat[l,4] <- as.character(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mut.rate[l])
                        mat[l,5] <- as.character(1/3)
                    }
                rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@locus.params <- as.data.frame(mat)
            }
})
