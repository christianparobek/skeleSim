#
# expressions to make a skelesim object from input$xxxx variables
#
# in general, wait till one of the inputs is changed and then alter the ssClass component of rValues


####general skelesim parameters
####This section updates the ssClass reactive when inputs change

observeEvent(input$title, {
    rValues$ssClass@title <- input$title
})
observeEvent(input$date, {
    rValues$ssClass@date <- as.POSIXct(input$date) #do we need POSIX dates?
})
observeEvent(input$quiet, {
    rValues$ssClass@quiet <- input$quiet
})
             
observeEvent(input$coalescent,{
    rValues$ssClass@simulator.type <- ifelse(input$coalescent,"c","f")
    rValues$ssClass@simulator <- ifelse(input$coalescent,"fsc","rmw")
    for (s in 1:length(rValues$ssClass@scenarios))
        if (rValues$ssClass@simulator.type=="c")
            {
                rValues$ssClass@scenarios[[s]]@simulator.params <-
                    fastsimcoalInit()
            } else {
                rValues$ssClass@scenarios[[s]]@simulator.params <-
                    rmetasimInit()
                rValues$ssClass@scenarios[[s]]@simulator.params@num.alleles <- rep(1,rValues$ssClass@scenarios[[s]]@num.loci)
                rValues$ssClass@scenarios[[s]]@simulator.params@allele.freqs <- vector("list",rValues$ssClass@scenarios[[s]]@num.loci)
                rValues$ssClass@scenarios[[s]]@simulator.params@num.gen <- 50
            }
})

observeEvent(input$reps, {
    rValues$ssClass@num.reps <- input$reps
})
observeEvent(input$wd, {
    rValues$ssClass@wd <- input$wd
                    })


###### scenario parameter updates
###### more complex, because it requires tracking which scenario
######

observeEvent(input$scenarioNumber,
             {
                 rValues$scenarioNumber <- input$scenarioNumber
             })

#observeEvent(input$migrationNumber,
#             {
#                 rValues$migrationNumber <- input$migrationNumber
#             })


observeEvent(input$numpopsTxt,
             {
                 numpop <- suppressWarnings(as.numeric(input$numpopsTxt))
                 if (!is.na(numpop)) 
                     rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops <- numpop
                 mig.mat()
             })
observeEvent(input$numloci,
             {
                 rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci <- input$numloci
                 #### rmetasim addition
                 if (rValues$ssClass@simulator.type=="f")
                 {
                     print("resetting the number of alleles per locus vector")
                     navec <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.alleles
                     print(navec)
                     #could test navec and append or shrink.  right now, we hose all num alleles values if the number of loci changes
                     if (length(navec)!=input$numloci) 
                         navec <- rep(1,input$numloci)
                     print(navec)
                     navec[is.na(navec)] <- 1
                     rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.alleles <- navec
                 }
             })

observeEvent(input$loctype,
             {
                 rValues$ssClass@scenarios[[rValues$scenarioNumber]]@locus.type <- input$loctype
             })

observeEvent(input$seqlen,
             {
                 rValues$ssClass@scenarios[[rValues$scenarioNumber]]@sequence.length <- input$seqlen
             })

observeEvent(input$migModel,
             {
                 rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migModel <- input$migModel
                 mig.mat()
             })

observeEvent(input$migRate,
             {
                 rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migRate <- input$migRate
                 mig.mat()
             })

observeEvent(input$rows,
             {
                 rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$rows <- input$rows
                 mig.mat()
             })

observeEvent(input$cols,
             {
                 rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$cols <- input$cols
                 mig.mat()
             })
observeEvent(input$distfun,
             {
                 rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$distfun <- input$distfun
                 mig.mat()
             })

####simcoal parameters
observeEvent(input$specScenNumber,
             {
                 rValues$scenarioNumber <- input$specScenNumber
             })

observeEvent(input$infSiteModel,
             {
                 if (rValues$ssClass@simulator.type=="c")
                     rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@inf.site.model <- input$infSiteModel
             })

observeEvent(input$fscexec,
             {
                 if (rValues$ssClass@simulator.type=="c")
                     rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@fastsimcoal.exec <- input$fscexec
             })

observeEvent(input$stvec,
             {
                 if (rValues$ssClass@simulator.type=="c")
                     rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@sample.times <- c(input$stvec)
             })

observeEvent(input$grvec,
             {
                 if (rValues$ssClass@simulator.type=="c")
                     rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@growth.rate <- c(input$grvec)
             })

#######rmetasim parameters
observeEvent(input$stages,
{
    if (rValues$ssClass@simulator.type=="f")
    {
        rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.stgs <- input$stages
        if ((is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@surv.matr)) || (dim(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@surv.matr)[1]!=input$stages))
        {
            rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@surv.matr <- matrix(0,input$stages,input$stages)
            rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@repr.matr <- matrix(0,input$stages,input$stages)
            rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@male.matr <- matrix(0,input$stages,input$stages)
        }
    
    }
})

observeEvent(input$self,
{
    if (rValues$ssClass@simulator.type=="f")
    {
         rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@selfing <- input$self
    }
})




##############################################################################################

#### this section updates the input boxes if a rValues$ssClass updates 
####
#######this observer is intended to run any time ssClass changes
####### it needs to be able to update inputs

observeEvent(rValues$ssClass,{
    if (!is.null(rValues$ssClass@title))
        updateTextInput(session,"title",value=rValues$ssClass@title)
    if (!is.null(rValues$ssClass@date))
        updateDateInput(session,"date",value=rValues$ssClass@date)
    if (!is.null(rValues$ssClass@quiet))
        updateCheckboxInput(session,"quiet",value=rValues$ssClass@quiet)
    if (!is.null(rValues$ssClass@simulator.type)){##sets a bunch of downstream parameters based on simulation type
        updateCheckboxInput(session,"coalescent",value=ifelse(rValues$ssClass@simulator.type=="c",T,F))
        output$simulator <- renderText({paste("Simulator:",rValues$ssClass@simulator)})
        if (rValues$ssClass@simulator.type=="c") rValues$ssClass@sim.func <- fsc.run else rValues$ssClass@sim.func <- rms.run
        output$simfunc <- renderText({paste("Simulator function:",ifelse(rValues$ssClass@simulator.type=="c","fsc.run","rms.run"))})
        if (rValues$ssClass@simulator.type=="c") rValues$ssClass@sim.check.func <- fsc.scenarioCheck else rValues$ssClass@sim.check.func <- rms.scenarioCheck
        
    }
    if (!is.null(rValues$ssClass@num.reps))
        updateNumericInput(session,"reps",value=rValues$ssClass@num.reps)
    if (!is.null(rValues$ssClass@wd))
        {
            if (is.null(supportValues$simroot)) {supportValues$simroot <- "."}
            output$simpath <- renderText({
                paste("Complete path for simulations to be executed:",paste(supportValues$simroot,rValues$ssClass@wd,sep="/"))
                  })
        }
##scenarios #respect the scenarioNumber!
    
    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops))
        updateTextInput(session,"numpopsTxt",value=paste(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops))

    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci))
        updateNumericInput(session,"numloci",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci)

    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@locus.type))
        updateNumericInput(session,"loctype",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@locus.type)

    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@sequence.length))
        updateNumericInput(session,"seqlen",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@sequence.length)

    ### needed for keeping track of how matrices are built in different scenarios
    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migModel))
        updateSelectInput(session,"migModel",selected=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migModel)
    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migRate))
        updateNumericInput(session,"migRate",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migRate)
    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$rows))
        updateNumericInput(session,"rows",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$rows)
    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$cols))
        updateNumericInput(session,"cols",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$cols)
    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$distfun))
        updateTextInput(session,"distfun",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$distfun)

####  this is the fastsimcoal updater
    if (input$coalescent)
        {
            if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@inf.site.model))
                {
                    updateCheckboxInput(session,"infSiteModel",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@inf.site.model)
                }
        } else { ########this is for rmetasim
#            if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@selfing))
#                updateNumericInput(session,"self",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@selfing)
            if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.gen))
                updateNumericInput(session,"gens",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.gen)
            if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.stgs))
                updateNumericInput(session,"stages",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.stgs)
        }
})

###change stuff if the scenario number changes
### this is a central 'function' that has grown organically.  In other words, its a mess.
###
observeEvent(rValues$scenarioNumber,
             {
                 if (!scenario.exists()) 
                     {
                         rValues$ssClass@scenarios <- c(rValues$ssClass@scenarios,rValues$ssClass@scenarios[1])
                     }                 

################ migration
###################
                 
                 if (input$coalescent) #simcoal
                     {
                         if (!is.null(rValues$ssClass@scenarios[[rValues$lstScenario]]@simulator.params))
                             rValues$ssClass@scenarios[[rValues$lstScenario]]@simulator.params@hist.ev <- rValues$history

                         if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@hist.ev))
                             {
#                                 print ("maybe should rewrite history?")
                                rValues$history <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@hist.ev
                             } else {
                                   rValues$history <- NULL
                             }
                     }


                 rValues$lstScenario <- rValues$scenarioNumber

                 ######## update the scenario input boxes
                 updateNumericInput(session,"scenarioNumber",value=rValues$scenarioNumber)
                 updateNumericInput(session,"specScenNumber",value=rValues$scenarioNumber)
             },priority=100)



###this observeEvent makes sure that the migration matrix stored in the reactive class is updated continuously
###very important!!!
observeEvent(input$migmat,{
    rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration[[1]] <- input$migmat
})

observeEvent(input$psvec,{
    rValues$ssClass@scenarios[[rValues$scenarioNumber]]@pop.size <- c(input$psvec)
})

observeEvent(input$ssvec,{
    rValues$ssClass@scenarios[[rValues$scenarioNumber]]@sample.size <- c(input$ssvec)
})

observeEvent(input$mutvec,{
    rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mut.rate <- c(input$mutvec)
})


### simcoal history updating
observeEvent(hst(),{
#    print("hst() observEvent")
    if (rValues$ssClass@simulator.type=="c")
        rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@hist.ev <- as.matrix(hst())
})
