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
    rValues$ssClass@simulator <- ifelse(input$coalescent,"fsc","rms")
    for (s in 1:length(rValues$ssClass@scenarios))
        if (rValues$ssClass@simulator.type=="c")
            {
                rValues$ssClass@scenarios[[s]]@simulator.params <-
                    fastsimcoalInit(rValues$ssClass@scenarios[[s]]@num.pops)
            } else {
                rValues$ssClass@scenarios[[s]]@simulator.params <-
                    rmetasimInit(rValues$ssClass@scenarios[[s]]@num.pops)
                rValues$ssClass@scenarios[[s]]@simulator.params@num.alleles <- rep(1,rValues$ssClass@scenarios[[s]]@num.loci)
                rValues$ssClass@scenarios[[s]]@simulator.params@allele.freqs <- vector("list",rValues$ssClass@scenarios[[s]]@num.loci)
                rValues$ssClass@scenarios[[s]]@simulator.params@num.gen <- 50
            }
})

observeEvent(input$reps, {
    rValues$ssClass@num.sim.reps <- as.numeric(floor(input$reps))
    updateNumericInput(session,"reps",value = rValues$ssClass@num.sim.reps)
})

observeEvent(input$NumPermReps,{
    rValues$ssClass@num.perm.reps <- floor(input$NumPermReps)
    updateNumericInput(session,"NumPermReps",value = rValues$ssClass@num.perm.reps)
})

observeEvent(input$analysesReq, {
    vec <- input$analysesReq
    if (length(vec)>0)
    {
        reqvec <- c("Global"=F,"Pairwise"=F,"Locus"=F)
        reqvec[which(names(reqvec)%in%input$analysesReq)] <- T
        rValues$ssClass@analyses.requested <- reqvec
    } else {
        rValues$ssClass@analyses.requested<- c("Global"=F,"Pairwise"=F,"Locus"=F)
    }


})

observeEvent(input$wd, {
    rValues$ssClass@wd <- input$wd
                    })


###### scenario parameter updates
###### more complex, because it requires tracking which scenario
######

observeEvent(input$scenarioNumber,
             {
                 if (req(input$scenarioNumber)>0)
                 {
                     rValues$scenarioNumber <- floor(input$scenarioNumber)
                 }
                 updateNumericInput(session,"scenarioNumber",value=rValues$scenarioNumber)
             })

#observeEvent(input$migrationNumber,
#             {
#                 rValues$migrationNumber <- input$migrationNumber
#             })


observeEvent(input$numpopsTxt,
             {
                 numpop <- suppressWarnings(floor(as.numeric(input$numpopsTxt)))
                 updateTextInput(session,"numpopsTxt",value=numpop)
                 if (!is.na(numpop))
                     rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops <- numpop
                 if (req(rValues$ssClass@simulator.type)=="c")
                     {
                         rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@growth.rate <- rep(0,rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops)
                             
                         rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@sample.times <- as.integer(floor(rep(0,rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops)))
                             
                     }
                 rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration[[rValues$migrationNumber +1]] <- mig.mat()
             })

observeEvent(input$numloci,
             {
                 if (!is.na(input$numloci))
                     if (input$numloci>0)
                     {
                         rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci <- floor(input$numloci)
                     }
                 updateNumericInput(session,"numloci",
                                    value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci)
                 #### rmetasim addition
                 if (rValues$ssClass@simulator.type=="f")
                 {
                     if (debug()) print("resetting the number of alleles per locus vector")
                     navec <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.alleles
#                     print(navec)
                     #could test navec and append or shrink.  right now, we hose all num alleles values if the number of loci changes
                     if (length(navec)!=floor(input$numloci))
                         navec <- rep(1,floor(input$numloci))
#                     print(navec)
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
                 rValues$ssClass@scenarios[[rValues$scenarioNumber]]@sequence.length <- floor(input$seqlen)
                 updateNumericInput(session,"seqlen",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@sequence.length)
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
                 rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$rows <- floor(input$rows)
                 updateNumericInput(session,"rows",value= rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$rows)
                 mig.mat()
             })

observeEvent(input$cols,
             {
                 rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$cols <- input$cols
                 updateNumericInput(session,"cols",value= rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$cols)
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
                 if (!is.na(input$specScenNumber))
                     if (input$specScenNumber>0)
                     {
                         rValues$scenarioNumber <- as.integer(floor(input$specScenNumber))
                     }
                 updateNumericInput(session,"specScenNumber",value=rValues$scenarioNumber)

             })

observeEvent(input$fscexec,
             {
                 if (rValues$ssClass@simulator.type=="c")
                 {
                     if (debug()) print(input$fscexec)
                     if (input$fscexec!="")
                         rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@fastsimcoal.exec <- input$fscexec
                 }
             })

observeEvent(input$stvec,
             {
                 if (rValues$ssClass@simulator.type=="c")
                     rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@sample.times <- c(as.integer(floor(input$stvec)))
             })

observeEvent(input$grvec,
             {
                 if (rValues$ssClass@simulator.type=="c")
                     rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@growth.rate <-
                     c(input$grvec)
             })

#######rmetasim parameters
observeEvent(input$stages,
{
    if (rValues$ssClass@simulator.type=="f")
    {
        
        rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.stgs <- floor(input$stages)
        updateNumericInput(session,"stages",
                           value= rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.stgs)
        if ((is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@surv.matr)) || (dim(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@surv.matr)[1]!=floor(input$stages)))
        {
            rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@surv.matr <- matrix(0,floor(input$stages),floor(input$stages))
            rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@repr.matr <- matrix(0,floor(input$stages),floor(input$stages))
            rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@male.matr <- matrix(0,floor(input$stages),floor(input$stages))
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
    if (!is.null(rValues$ssClass@num.perm.reps))
        updateNumericInput(session,"NumPermReps",value=rValues$ssClass@num.perm.reps)
    if (!is.null(rValues$ssClass@quiet))
        updateCheckboxInput(session,"quiet",value=rValues$ssClass@quiet)
    if (!is.null(rValues$ssClass@simulator.type)){##sets a bunch of downstream parameters based on simulation type
        updateCheckboxInput(session,"coalescent",value=ifelse(rValues$ssClass@simulator.type=="c",T,F))
        output$simulator <- renderText({paste("Simulator:",rValues$ssClass@simulator)})
        if (rValues$ssClass@simulator.type=="c") rValues$ssClass@sim.func <- fsc.run else rValues$ssClass@sim.func <- rms.run
        output$simfunc <- renderText({paste("Simulator function:",ifelse(rValues$ssClass@simulator.type=="c","fsc.run","rms.run"))})
        if (rValues$ssClass@simulator.type=="c") rValues$ssClass@sim.check.func <- fsc.scenarioCheck else rValues$ssClass@sim.check.func <- rms.scenarioCheck

    }
    if (!is.null(rValues$ssClass@num.sim.reps))
        updateNumericInput(session,"reps",value=rValues$ssClass@num.sim.reps)

    if (!is.null(rValues$ssClass@analyses.requested))
        updateCheckboxGroupInput(session,"analysesReq",
                                 names(rValues$ssClass@analyses.requested)[rValues$ssClass@analyses.requested])


        if (!is.null(rValues$ssClass@wd))
        {
            if (is.null(supportValues$simroot)) {supportValues$simroot <- "."}
            output$simpath <- renderText({
                paste("Complete path for simulation 'scratch' directory:",paste(supportValues$simroot,rValues$ssClass@wd,sep="/"))
                  })
        }
##scenarios #respect the scenarioNumber!

    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops))
    {
        updateTextInput(session,"numpopsTxt",
                        value=paste(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops))
        rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration[[rValues$migrationNumber+1]] <- mig.mat()
    }

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

        } else { ########this is for rmetasim

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
                         rValues$ssClass@scenarios <-
                             c(rValues$ssClass@scenarios,rValues$ssClass@scenarios[1])
                     }

################ migration
###################

                 if (input$coalescent) #simcoal
                     {
#                         if (!is.null(rValues$ssClass@scenarios[[rValues$lstScenario]]@simulator.params))
#                             rValues$ssClass@scenarios[[rValues$lstScenario]]@simulator.params@hist.ev <- rValues$history

                         if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@hist.ev))
                             {
#                                 print ("maybe should rewrite history?")
#                                rValues$history <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@hist.ev
                             } else {
#                                   rValues$history <- NULL
                             }
                     } else { #rmetasim stuff
                         
                     }


                 rValues$lstScenario <- rValues$scenarioNumber

                 ######## update the scenario input boxes
                 updateNumericInput(session,"scenarioNumber",value=rValues$scenarioNumber)
                 updateNumericInput(session,"specScenNumber",value=rValues$scenarioNumber)
             },priority=100)



###this observeEvent makes sure that the migration matrix stored in the reactive class is updated continuously
###very important!!!
observeEvent(input$migmat,{
    rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration[[rValues$migrationNumber+1]] <- input$migmat
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

#this is a change in the migration matrix number
#changing will require 1) making sure the new one exists or 2) create it, and 3) set variables
observeEvent(input$mignum,{
    mn <- floor(input$mignum)
    if (mn!=rValues$migrationNumber)
    {
        migs <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration
        if (length(migs)<(mn+1)) #need another x matrices
        {
            migs <- c(migs,rep(migs[length(migs)],(mn+1)-length(migs))) #create new ones out of the last matrix
        }
        rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration <- migs
        rValues$migrationNumber <- mn
        updateNumericInput(session,"mignum",value=mn)
        if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migModel))
        {
            rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migModel <- "user"
            updateSelectInput(session,"migModel",selected="user")
        }

    }
})

output$numMigMats <- renderText({
    req(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration)
    if (length(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration)==1) mtxt <- "matrix" else mtxt <- "matrices"
    paste(length(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration),"migration",mtxt,"defined currently (indexed from '0')")
    })


### Simcoal history updating
observeEvent(hst(),{
    if (req(rValues$ssClass@simulator.type)=="c")
    {
        if (debug()) print("about to assign hist.ev.  Current value:")
        if (debug()) print(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@hist.ev)
        rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@hist.ev <- as.matrix(hst())
    }
})

