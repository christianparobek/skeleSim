#
# expressions to create a skelesim object from input$xxxx variables 
#
# maintain and build a list of scenarios
# list of scenarios is appended if the scenario number
# is incremented

scenario.observe <- observe({     
    sn <- input$scenarioNumber

    if (is.null(sn)) {sn <- 1}   #at least one scenario
    print(paste("sn",sn))
    locscen <- scenarios

    print(str(locscen))
    
    print(paste("len scen",length(scenarios)))
    print(paste("len locscen",length(locscen)))
    if (!is.null(input$numpops))
        {
    ####here is the situation where there is a new scenario added
    if (sn<=length(locscen))
        {
            print("enter cond true")
            locscen[[sn]] <- new("scenario.params")
            locscen[[sn]]@num.pops <- input$numpops
            locscen[[sn]]@migration <- inmat()
                        print("so far")

            print("before rep")
            locscen[[sn]]@pop.size <- rep(100,input$numpops) #hard-coded need to fix
            locscen[[sn]]@sample.size <- rep(25,input$numpops) #hard-coded need to fix
            locscen[[sn]]@locus.type <- "sequence"            #hard-coded
            locscen[[sn]]@sequence.length <- 400            #hard-coded
            print("after rep")
            locscen[[sn]]@num.loci <- input$numloci
            locscen[[sn]]@mut.rate <- input$mutRate
        } else  {
            print(paste("sn",sn))
            locscen <- c(locscen,locscen[rep(1,(length(locscen)-sn))])
            for (s in (length(locscen)+1):sn)
                {
                    locscen[[s]] <- locscen[[1]]
                }
        }
    scenarios<<-locscen
}
})

scenario.return <- reactive( {scenarios[[input$scenarioNumber]]} )

output$scenDebug <- renderText({scenario.return()})


skelesim.class <- reactive({
    
params <- new("skeleSim.params")
params@title <- input$title
params@date <- input$date
params@quiet <- input$quiet
params@question <- "n"      #not sure what to do here yet
params@simulator.type <- ifelse(input$coalescent,"c","i")
params@simulator <- "fsc"   #this is hard-coded but might need some work
params@num.reps <- input$reps
params@timing <- input$timing
params@sim.func <- fsc.run  #hard-coded right now
params@wd <- input$wd
})

