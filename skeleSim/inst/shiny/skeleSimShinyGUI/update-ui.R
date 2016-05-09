##############################################################################################

#### this section updates the input boxes if a rValues$ssClass updates
####
#######this observer is intended to run any time ssClass changes
####### it needs to be able to update inputs

#observeEvent(rValues$ssClass,{

observeEvent(rValues$ssClass@title,{
    if (!is.null(rValues$ssClass@title))
        updateTextInput(session,"title",value=rValues$ssClass@title)
})
observeEvent(rValues$ssClass@date,{
    if (!is.null(rValues$ssClass@date))
        updateDateInput(session,"date",value=rValues$ssClass@date)
})
observeEvent(rValues$ssClass@num.perm.reps,{
    if (!is.null(rValues$ssClass@num.perm.reps))
        updateNumericInput(session,"NumPermReps",value=rValues$ssClass@num.perm.reps)
})
observeEvent(rValues$ssClass@quiet,{
    if (!is.null(rValues$ssClass@quiet))
        updateCheckboxInput(session,"quiet",value=rValues$ssClass@quiet)
})
observeEvent(rValues$ssClass@simulator.type,{
    if (!is.null(rValues$ssClass@simulator.type)){##sets a bunch of downstream parameters based on simulation type
        updateCheckboxInput(session,"coalescent",value=ifelse(rValues$ssClass@simulator.type=="c",T,F))
        output$simulator <- renderText({paste("Simulator:",rValues$ssClass@simulator)})
        if (rValues$ssClass@simulator.type=="c") rValues$ssClass@sim.func <- fsc.run else rValues$ssClass@sim.func <- rms.run
        output$simfunc <- renderText({paste("Simulator function:",ifelse(rValues$ssClass@simulator.type=="c","fsc.run","rms.run"))})
        if (rValues$ssClass@simulator.type=="c") rValues$ssClass@sim.check.func <- fsc.scenarioCheck else rValues$ssClass@sim.check.func <- rms.scenarioCheck
        
    }
})

observeEvent(rValues$ssClass@num.sim.reps,{   
    if (!is.null(rValues$ssClass@num.sim.reps))
        updateNumericInput(session,"reps",value=rValues$ssClass@num.sim.reps)
})

observeEvent(rValues$ssClass@analyses.requested,{
    if (!is.null(rValues$ssClass@analyses.requested))
        updateCheckboxGroupInput(session,"analysesReq",
                                 names(rValues$ssClass@analyses.requested)[rValues$ssClass@analyses.requested])
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
    
#observeEvent(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci,{
#    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci))
#        updateNumericInput(session,"numloci",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci)
#})

#observeEvent(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@locus.type,{
#    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@locus.type))
#        updateNumericInput(session,"loctype",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@locus.type)
#})

#observeEvent(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@sequence.length,{
#    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@sequence.length))
#        updateNumericInput(session,"seqlen",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@sequence.length)
#})
### needed for keeping track of how matrices are built in different scenarios
#observeEvent(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migModel,{
#    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migModel))
#        updateSelectInput(session,"migModel",selected=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migModel)
#})

#observeEvent(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migRate,{
#    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migRate))
#        updateNumericInput(session,"migRate",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migRate)
#})

#observeEvent(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$rows,{
#    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$rows))
#        updateNumericInput(session,"rows",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$rows)
#})
#observeEvent(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$cols,{
#    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$cols))
#        updateNumericInput(session,"cols",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$cols)
#})

#observeEvent(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$distfun,{
#    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$distfun))
#        updateTextInput(session,"distfun",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$distfun)
#
#})
        
observeEvent(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params,{
    if (rValues$ssClass@simulator.type!="c") #rmetasim, for now
    {
        updateNumericInput(session,"gens",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.gen)
        updateNumericInput(session,"stages",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.stgs)
    }
})


updateUIs <- function()
{
                    if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops))
                     updateTextInput(session,"numpopsTxt",
                                     value=paste(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops))

                  if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci))
                      updateNumericInput(session,"numloci",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci)
                 if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@locus.type))
                     updateNumericInput(session,"loctype",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@locus.type)
                 if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@sequence.length))
                     updateNumericInput(session,"seqlen",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@sequence.length)
                 if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migModel))
                     updateSelectInput(session,"migModel",selected=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migModel)
                 if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$rows))
                     updateNumericInput(session,"rows",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$rows)
                 if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$cols))
                         updateNumericInput(session,"cols",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$cols)
                 if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$distfun))
                     updateTextInput(session,"distfun",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$distfun)
                 if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migRate))
                     updateNumericInput(session,"migRate",value=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migRate)
}
