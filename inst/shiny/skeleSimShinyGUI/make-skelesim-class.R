#
# expressions to make a skelesim object from input$xxxx variables
#
# in general, wait till one of the inputs is changed and then alter the ssClass component of rValues


####general skelesim parameters
####This section updates the ssClass reactive when inputs change

observeEvent(input$title, {
    rValues$ssClass@title <- isolate(input$title)
})
observeEvent(input$date, {
    rValues$ssClass@date <- isolate(as.POSIXct(input$date)) #do we need POSIX dates?
})
observeEvent(input$quiet, {
    rValues$ssClass@quiet <- isolate(input$quiet)
})

observeEvent(input$coalescent,{
    rValues$ssClass@simulator.type <- ifelse(isolate(input$coalescent),"c","f")
    rValues$ssClass@simulator <- ifelse(isolate(input$coalescent),"fsc","rms")
    for (s in 1:length(rValues$ssClass@scenarios))
        if (rValues$ssClass@simulator.type=="c")
            {
                if (is.null(rValues$ssClass@scenarios[[s]]@simulator.params)|class(rValues$ssClass@scenarios[[s]]@simulator.params)!="fastsimcoal.params")
                {
                    rValues$ssClass@scenarios[[s]]@simulator.params <-
                        isolate(fastsimcoalInit(rValues$ssClass@scenarios[[s]]@num.pops))
                }
            } else {
                if (is.null(rValues$ssClass@scenarios[[s]]@simulator.params)|class(rValues$ssClass@scenarios[[s]]@simulator.params)!="rmetasim.params")
                {
                    rValues$ssClass@scenarios[[s]]@simulator.params <-
                        rmetasimInit(rValues$ssClass@scenarios[[s]]@num.pops)
                    rValues$ssClass@scenarios[[s]]@simulator.params@num.alleles <- isolate(rep(1,rValues$ssClass@scenarios[[s]]@num.loci))
                    rValues$ssClass@scenarios[[s]]@simulator.params@allele.freqs <- isolate(vector("list",rValues$ssClass@scenarios[[s]]@num.loci))
                    rValues$ssClass@scenarios[[s]]@simulator.params@num.gen <- 50
                }
            }
})

observeEvent(input$reps, {
    rValues$ssClass@num.sim.reps <- isolate(as.numeric(floor(input$reps)))

})

observeEvent(input$NumPermReps,{
    rValues$ssClass@num.perm.reps <- isolate(floor(input$NumPermReps))

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
                     rValues$scenarioNumber <- isolate(floor(as.numeric(input$scenarioNumber)))
                 }
#                 updateNumericInput(session,"scenarioNumber",value=isolate(rValues$scenarioNumber))
             })



observeEvent(input$numpopsTxt,
             {
                 numpop <- suppressWarnings(floor(as.numeric(input$numpopsTxt)))
 #                updateTextInput(session,"numpopsTxt",value=numpop)
                 if (!is.na(numpop))
                     rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops <- numpop
                 if (req(rValues$ssClass@simulator.type)=="c")
                     {
                         rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@growth.rate <- rep(0,rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops)
                             
                         rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@sample.times <- as.integer(floor(rep(0,rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops)))
                             
                     }
                 rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration[[rValues$migrationNumber +1]] <- mig.mat()
             },priority=-1)



observeEvent(input$numloci,
             {
                 numloci <- as.numeric(input$numloci)
                 if (!is.na(numloci))
                 {
                     if (numloci>0)
                     {
                         rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci <- floor(numloci)
                         if (length(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mut.rate)
                             !=
                             rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci)
                         {
                             rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mut.rate <- rep(1e-5,
                                                                                                 rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci)
                         }
                     }
                 }
                 #### rmetasim addition
                 if (!is.na(numloci))
                     if (rValues$ssClass@simulator.type=="f")
                     {
                         if (debug()) print("resetting the number of alleles per locus vector")
                         navec <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.alleles
####                     print(navec)
####could test navec and append or shrink.  right now, we hose all num alleles values if the number of loci changes
                         if (length(navec)!=floor(numloci))
                             navec <- rep(1,floor(numloci))
###                     print(navec)
                         navec[is.na(navec)] <- 1
                         rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.alleles <- navec
                     }
             },priority=-1)

observe({
    if (rValues$ssClass@simulator.type=="f")
        output$focalLoc <- renderUI({
            numericInput("focalLoc","Adjust allele frequencies for locus",value=1,min=1,max=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci)
        })
})

observeEvent(input$loctype,
             {
                 req(input$loctype)
                 if (input$loctype=="sequence")
                 {
                     rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci <- 1
                     updateTextInput(session,"numloci",value="1")
                 } else {

                 }
                 
                 rValues$ssClass@scenarios[[rValues$scenarioNumber]]@locus.type <- input$loctype
             })

observeEvent(input$seqlen,
             {
                 rValues$ssClass@scenarios[[rValues$scenarioNumber]]@sequence.length <- floor(input$seqlen)
                 #updateNumericInput(session,"seqlen",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@sequence.length)
             })

observeEvent(input$migModel,
             {
                 rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migModel <- input$migModel
                # mig.mat()
             })

observeEvent(input$migRate,
             {
                 rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migRate <- input$migRate
                 #mig.mat()
             })

observeEvent(input$rows,
             {
                 rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$rows <- floor(input$rows)
#                 updateNumericInput(session,"rows",value= rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$rows)
                 mig.mat()
             })

observeEvent(input$cols,
             {
                 rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$cols <- input$cols
#                 updateNumericInput(session,"cols",value= rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$cols)
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
                         rValues$scenarioNumber <- isolate(as.integer(floor(as.numeric(input$specScenNumber))))
                     }
 #                updateNumericInput(session,"specScenNumber",value=isolate(rValues$scenarioNumber))

             })

observeEvent(input$fscexec,
             {
                 if (rValues$ssClass@simulator.type=="c")
                 {
                     if (debug()) print(input$fscexec)
                     if (input$fscexec!="")
                         rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@fastsimcoal.exec <- isolate(input$fscexec)
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
        rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.stgs <- isolate(floor(input$stages))
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

observeEvent(input$gens,{
    if (rValues$ssClass@simulator.type=="f")
    {
        if (req(input$gens)>=1)
            rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.gen <- input$gens
    }
    
})



##############################################################################################

#### this section updates the input boxes if a rValues$ssClass updates
####
#######this observer is intended to run any time ssClass changes
####### it needs to be able to update inputs

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

                 sanityChecks()

                 updateUIs()
 
             },priority=100)



###this observeEvent makes sure that the migration matrix stored in the reactive class is updated continuously
###very important!!!
observeEvent(input$migmat,{
    rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration[[rValues$migrationNumber+1]] <- input$migmat
})

observeEvent(input$psvec,{
if (debug())    print("in observevent psvec")
if (debug())      print(rValues$scenarioNumber)
if (debug())      print(input$psvec)
if (debug())      print(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@pop.size)
if (debug())      print(length(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@pop.size))
    if (length(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@pop.size) != rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops)
    {
        rValues$ssClass@scenarios[[rValues$scenarioNumber]]@pop.size <-
            rep(100,rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops)
    } else {
        rValues$ssClass@scenarios[[rValues$scenarioNumber]]@pop.size <- c(input$psvec)
    }
    
if (debug())      print(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@pop.size)
if (debug())      print("in observevent psvec")
})

observeEvent(input$ssvec,{
    rValues$ssClass@scenarios[[rValues$scenarioNumber]]@sample.size <- c(input$ssvec)
})

observeEvent(input$mutvec,{

    req(input$mutvec)

if (debug())      print("in mutvec observeEvent")
if (debug())      print(rValues$scenarioNumber)
if (debug())      print(input$mutvec)
if (debug())      print(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mut.rate)
   diff <- length(input$mutvec)-length(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mut.rate)
if (debug())      print(diff)

#    if (diff==0) #lengths are good. replace
    
        rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mut.rate <- c(input$mutvec)
        
#    else if (diff<0) #input short
#        rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mut.rate <- c(input$mutvec,rep(0.0001,abs(diff)))
#    else #input long
#        {
#            rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mut.rate <-
#                c(input$mutvec)[1:(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci)]
#        }
if (debug())      print(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mut.rate)
if (debug())      print("leaving observerevent mutvec")

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
        updateNumericInput(session,"mignum",value=isolate(mn))
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


###scenario adding
observeEvent(input$addScenario,{
    rValues$scenarioNumber <- isolate(length(rValues$ssClass@scenarios)+1)
})
