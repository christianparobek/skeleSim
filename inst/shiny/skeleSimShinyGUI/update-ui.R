##############################################################################################

#### this section updates the input boxes if a rValues$ssClass updates
####
#######this observer is intended to run any time ssClass changes
####### it needs to be able to update inputs

#observeEvent(rValues$ssClass,{

observeEvent(rValues$ssClass@title,{
    if (!is.null(rValues$ssClass@title))
        updateTextInput(session,"title",value=isolate(rValues$ssClass@title))
})
observeEvent(rValues$ssClass@date,{
    if (!is.null(rValues$ssClass@date))
        updateDateInput(session,"date",value=isolate(rValues$ssClass@date))
})
observeEvent(rValues$ssClass@num.perm.reps,{
    if (!is.null(rValues$ssClass@num.perm.reps))
        updateNumericInput(session,"NumPermReps",value=isolate(rValues$ssClass@num.perm.reps))
})
observeEvent(rValues$ssClass@quiet,{
    if (!is.null(rValues$ssClass@quiet))
        updateCheckboxInput(session,"quiet",value=isolate(rValues$ssClass@quiet))
})
observeEvent(rValues$ssClass@simulator.type,{
    if (!is.null(rValues$ssClass@simulator.type)){##sets a bunch of downstream parameters based on simulation type
        updateCheckboxInput(session,"coalescent",value=isolate(ifelse(rValues$ssClass@simulator.type=="c",T,F)))
        output$simulator <- renderText({paste("Simulator:",rValues$ssClass@simulator)})
        if (rValues$ssClass@simulator.type=="c") rValues$ssClass@sim.func <- fsc.run else rValues$ssClass@sim.func <- rms.run
        output$simfunc <- renderText({paste("Simulator function:",ifelse(rValues$ssClass@simulator.type=="c","fsc.run","rms.run"))})
        if (rValues$ssClass@simulator.type=="c") rValues$ssClass@sim.check.func <- fsc.scenarioCheck else rValues$ssClass@sim.check.func <- rms.scenarioCheck
        
    }
})

observeEvent(rValues$ssClass@num.sim.reps,{   
    if (!is.null(rValues$ssClass@num.sim.reps))
        updateNumericInput(session,"reps",value=isolate(rValues$ssClass@num.sim.reps))
})

observeEvent(rValues$ssClass@analyses.requested,{
    if (!is.null(rValues$ssClass@analyses.requested))
        updateCheckboxGroupInput(session,"analysesReq",
                                 isolate(names(rValues$ssClass@analyses.requested)[rValues$ssClass@analyses.requested]))
})
    
observeEvent(rValues$ssClass@wd,{
    if (!is.null(rValues$ssClass@wd))
    {
        if (is.null(supportValues$simroot)) {supportValues$simroot <- "."}
        output$simpath <- renderText({
            paste("Complete path for simulation 'scratch' directory:",paste(supportValues$simroot,rValues$ssClass@wd,sep="/"))
        })
    }
})
    
    ##scenarios #respect the scenarioNumber!
    
observeEvent(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops,{
    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops))
    {
#        updateTextInput(session,"numpopsTxt",
#                        value=paste(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops))
        rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration[[rValues$migrationNumber+1]] <- mig.mat()
    }
})
    
        
observeEvent(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params,{
    if (rValues$ssClass@simulator.type!="c") #rmetasim, for now
    {
        updateNumericInput(session,"gens",value=isolate(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.gen))
        updateNumericInput(session,"stages",value=isolate(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.stgs))
    }
})


#
# this function updates a bunch of user interface items based on the rValues$ssClass entries
#

updateUIs <- function()
{


    ###genpars
    updateNumericInput(session,"reps",value = isolate(rValues$ssClass@num.sim.reps))
    updateNumericInput(session,"NumPermReps",value = isolate(rValues$ssClass@num.perm.reps))
    
    
                 ######## update the scenario input boxes
    choices=1:length(rValues$ssClass@scenarios)
    updateSelectInput(session,"scenarioNumber",choices=choices,selected=isolate(as.character(rValues$scenarioNumber)))
    updateSelectInput(session,"specScenNumber",choices=choices,selected=isolate(as.character(rValues$scenarioNumber)))
    

    
    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops))
        updateTextInput(session,"numpopsTxt",
                        value=isolate(paste(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops)))
    
    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci))
        updateNumericInput(session,"numloci",value=isolate(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci))
    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@locus.type))
        updateSelectInput(session,"loctype",selected=isolate(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@locus.type))
    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@sequence.length))
        updateNumericInput(session,"seqlen",value=isolate(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@sequence.length))
    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migModel))
        updateSelectInput(session,"migModel",selected=isolate(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migModel))
    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$rows))
        updateNumericInput(session,"rows",value=isolate(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$rows))
    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$cols))
        updateNumericInput(session,"cols",value=isolate(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$cols))
    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$distfun))
        updateTextInput(session,"distfun",value=isolate(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$distfun))
    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migRate))
        updateNumericInput(session,"migRate",value=isolate(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migRate))

###supposed to reset the executables section for fastsimcoal
    if (!is.null(rValues$ssClass@simulator))
        if (rValues$ssClass@simulator=="fsc")
            if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params))
                    {
                        sim.exec <-
                            c(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@fastsimcoal.exec,
                                      supportValues$simexec)
                        sim.exec <- sim.exec[!is.null(sim.exec)]
                        sim.exec <- basename(Sys.which(sim.exec)) #gets the names of the executables actually on machine
                        sim.exec <- unique(sim.exec)
                        sim.exec <- sim.exec[nchar(sim.exec)>2]
                        if (debug()) print("rendering UI for simexec")
                        if (debug()) print(sim.exec)
                        updateSelectInput(session,"fscexec",selected=sim.exec[1],
                                          choices=sim.exec)
                    }

    ###rmetasim
    if (req(rValues$ssClass@simulator)=="rms")
        {
            updateNumericInput(session,"stages",
                               value= isolate(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.stgs))
            nl = rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci
            updateNumericInput(session,"focalLoc",min=1,max=nl)

        }


    
    ##analysis results
    if (!is.null(rValues$ssClass@analysis.results))
    {
        updateSelectizeInput(session,"gstatsel",choices=isolate(global.stats(rValues$ssClass)),selected=isolate(global.stats(rValues$ssClass)),server=TRUE)

        updateSelectizeInput(session,"lstatsel",choices=isolate(locus.stats(rValues$ssClass)),selected=isolate(locus.stats(rValues$ssClass)),server=TRUE)

        updateSelectizeInput(session,"pstatsel",choices=isolate(pairwise.stats(rValues$ssClass)), selected=isolate(pairwise.stats(rValues$ssClass)),server=TRUE)

        updateNumericInput(session,"vizScenario",max=viz.scenarios(isolate(rValues$ssClass)))
    }
                    
}
