#
# expressions to create a skelesim object from input$xxxx variables
#
# maintain and build a list of scenarios
# list of scenarios is appended if the scenario number
# is incremented

sim.specific.update <- reactive({
    if (!is.null(ssClass@scenarios[[input$scenarioNumber]]))
        {
            if (is.null(ssClass@scenarios[[input$scenarioNumber]]@simulator.params))
                if (input$coalescent) #fast simcoal
                    {
                        simparam <- new("fastsimcoal.params")
                        simparam@hist.ev <- hst()
                    } else {
                        simparam <- new("rmetasim.params")
                    }
            simparam
        }
})


scenario.update <- reactive({
    ##### test for nulls in the reactives
    anynull <- FALSE
#    browser()
    if (is.null(input$scenarioNumber)&(!anynull)) anynull <- TRUE
    if (is.null(input$numpops)&(!anynull)) anynull <- TRUE
    if (is.null(input$numloci)&(!anynull)) anynull <- TRUE
    if (is.null(input$mut.rate)&(!anynull)) anynull <- TRUE
    if (is.null(input$migModel)&(!anynull)) anynull <- TRUE
    if (is.null(input$migRate)&(!anynull)) anynull <- TRUE

#    print (paste("anynull",anynull))
    locscen <- NULL
    if (!anynull)
        {
            sn <- input$scenarioNumber
            print(paste("sn first",sn))
            if (is.null(ssClass@scenarios))
                {
                    ssClass@scenarios <- list(new("scenario.params"))
                }
            locscen <- ssClass@scenarios
            print(paste("len locscen",length(locscen)))
            
                                        #        print(input$numpops)
####here is the situation where there is a new scenario added
            if (sn<=length(locscen))
                {
                    print("scenario number greater than the total number of scenarios to this point")
                    locscen[[sn]]@num.pops <- input$numpops
                    locscen[[sn]]@migration <- list(ssClass.scenario.mig())
                    locscen[[sn]]@mig.helper <- list(mig.model = input$migModel,
                                                     mig.rate  = input$migRate,
                                                     landrows  = input$rows,
                                                     landcols  = input$cols,
                                                     distfun   = input$distfun)
                    locscen[[sn]]@pop.size <- rep(100,input$numpops) #hard-coded need to fix
                    locscen[[sn]]@sample.size <- rep(25,input$numpops) #hard-coded need to fix
                    locscen[[sn]]@locus.type <- "sequence"            #hard-coded
                    locscen[[sn]]@sequence.length <- 400            #hard-coded
                    locscen[[sn]]@num.loci <- input$numloci
                    locscen[[sn]]@mut.rate <- input$mut.rate
                    locscen[[sn]]@simulator.params <- sim.specific.update()
                } else  {
                    print ("scenario number in the range of existing scenarios")
                    print(paste("scenario number in the scenario.observe observer:",sn))
                                        #                print(paste("length locscen2",length(locscen)))
                    oldlen <- length(locscen)
                    locscen <- c(locscen,locscen[rep(1,(sn-length(locscen)))])
                    for (s in (oldlen+1):sn)
                        {
                            locscen[[s]] <- locscen[[1]]
                        }
                }
                                        # browser(length(locscen)>0)
            print("scenarios:")
            print(str(locscen))
        }
    locscen
})



ssClass.update <- reactive({
    ssClass@title <<- input$title
    ssClass@date <<- as.POSIXct(Sys.time())
    ssClass@quiet <<- input$quiet
    ssClass@question <<- "n"      #not sure what to do here yet
    print(paste("coal",input$coalescent))
    if (!is.null(input$coalescent))
        {
            if (input$coalescent)
                {
                    ssClass@simulator.type <<- "c"
                    ssClass@sim.func <<- fsc.run
                } else {
                    ssClass@simulator.type <<- "f"
                    ssClass@sim.func <<- rms.run
                }
        }
    ssClass@simulator <<- "fsc"   #this is hard-coded but might need some work
    ssClass@num.reps <<- input$reps
    ssClass@wd <<- input$wd
    ssClass@scenarios <- scenario.update() 
                                        #print(ssClass)
    ssClass
})

observeEvent(ssClass.update(),{ssClass<<-ssClass.update()})
