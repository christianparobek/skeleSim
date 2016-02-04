#
# file to be included in server.R that specifies parts of the rmetasim backend
#



############
############ demography
############

observeEvent(callModule(matrixIn,"survmat",label="Survival Matrix",mat=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@surv.matr)(),{
    rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@surv.matr <- callModule(matrixIn,"survmat",label="Survival Matrix",mat=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@surv.matr)()
})

observeEvent(callModule(matrixIn,"repmat",label="Reproduction Matrix",mat=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@repr.matr)(),{
    rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@repr.matr <- callModule(matrixIn,"repmat",label="Reproduction Matrix",mat=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@repr.matr)()
})

observeEvent(callModule(matrixIn,"malemat",label="From col to row male contribution (defaults are fine)",mat=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@male.matr)(),{
    rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@male.matr <- callModule(matrixIn,"malemat",label="From col to row male contribution (defaults are fine)",mat=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@male.matr)()
})

########## 
output$lefkovitch <- renderTable({
    rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@surv.matr +
        rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@repr.matr
})

output$leading <- renderText({
    ev <- eigen(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@surv.matr +
        rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@repr.matr)$values
    paste("lambda",round(ev[1],2))
})


##in theory this observe event waits to see if the number of alleles at a locus changes and updates ssClass

numal.react <- reactive({
    vf=NULL
    if (req(rValues$ssClass@simulator.type=="f"))
    {
        nal <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.alleles
        if (!is.null(req(nal)))
        {
            nal[is.na(nal)] <- 1
            vf <- callModule(vectorIn,"numall",label="numalleles",
                             vec=nal)()
        }
    }
    vf
})
observeEvent(numal.react(),
                        {
                            rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.alleles <-
                                c(numal.react())
                        })


anums <- reactive({  #allele numbers pulled from reactive value.  this allows checks for existence
    if (req(rValues$ssClass@simulator.type)=="f")
        {
            rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.alleles
        } else {NULL}
})

observeEvent(anums(),
{
    print("running anums()")
    if (is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@allele.freqs))
    {
        rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@allele.freqs <- lapply(anums(),function(x){rep(1/x,x)})
    } else {
        if (length(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@allele.freqs)!=length(anums()))
        {
            print("setting allele freqs")
            rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@allele.freqs <- lapply(anums(),function(x){rep(1/x,x)})
        } else { #if there is a legitimate list and only changing the number of alleles for a locus be more smart
            aflst <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@allele.freqs #just to make easier to work with
            svec <- which(sapply(aflst,length) !=anums())
            for (i in svec) {aflst[[i]] <- rep(1/anums()[i],anums()[i])} #only replace elements whose allele numbers have changed
            rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@allele.freqs <- aflst
        }
    }
    
})


output$afreqs <- renderText({
    if (!is.null(req(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@allele.freqs)))
        paste(sapply(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@allele.freqs,function(x){paste(x,sep=", ")}),sep="\n")
       
})

#observe({
#    ANvecFn <- callModule(vectorIn,"numall",label="Allele num per loc",vec=
#})
