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
    mat <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@surv.matr +
        rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@repr.matr
    mat <- as.data.frame(mat)
    names(mat) <- paste0("From.",1:dim(mat)[2])
    rownames(mat) <- paste0("To.",1:dim(mat)[1])
    mat
})

output$leading <- renderText({
    ev <- eigen(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@surv.matr +
        rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@repr.matr)$values
    paste("discrete-time growth rate (lambda) implied by matrix above:",round(ev[1],2))
})


##in theory this observe event waits to see if the number of alleles at a locus changes and updates ssClass

#numal.react <- reactive({
#    vf=NULL
#    if (req(rValues$ssClass@simulator.type=="f"))
#    {
#        nal <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.alleles
#        if (!is.null(req(nal)))
#        {
#            nal[is.na(nal)] <- 1
#            vf <- callModule(vectorIn,"numall",label="Num. Alleles per locus (each element of vector corresponds to a locus)",
#                             vec=nal)()
#        }
#    }
#    vf
#})
#observeEvent(numal.react(),
#                        {
#                            rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.alleles <-
#                                c(numal.react())
#                        })


anums <- reactive({  #allele numbers pulled from reactive value.  this allows checks for existence
    if (req(rValues$ssClass@simulator.type)=="f")
        {
            rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.alleles
        } else {NULL}
})

num.alleles <- function() #a function to return a sensible number of alleles vector
    {
          if (debug()) print("running num.alleles()")
        if (req(rValues$ssClass@simulator.type=="f"))
        {
             nal <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.alleles
             if (!is.null(req(nal)))
             {
                 nal[is.na(nal)] <- 1
             }
             nal
        } else {c(1)}
    }

output$numall <- renderUI({
    if (debug()) print("creating allele num vector ui")

    matrixInput("numall","Number of alleles segregating at each locus at start of simulation",
                as.data.frame(t(num.alleles())))
})


observeEvent(input$numall,{
     rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.alleles <- c(input$numall)
})

observeEvent(anums(),
{
    if (debug()) print("running anums()")
    if ((req(rValues$ssClass@simulator.type)=="f"))
        if (is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@allele.freqs))
        {
            rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@allele.freqs <- lapply(anums(),function(x){rep(1/x,x)})
        } else   {
            if (length(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@allele.freqs)!=length(anums()))
            {
                if (debug()) print("setting allele freqs")
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


output$afreqLoc <- renderUI({
    if (debug())          print("creating allele freq  vector ui")
    if (!is.null(req(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@allele.freqs))&
        (!is.null(input$focalLoc)))
    {
        af <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@allele.freqs[[input$focalLoc]]
        matrixInput("afreqLoc",paste("Distribution of allele freqs at locus",input$focalLoc),
                    as.data.frame(matrix(af,nrow=1)))
    }
})

observeEvent(input$afreqLoc,{
    af = c(input$afreqLoc)
    af=abs(af)
    l = length(af)
    if (l>1)
        af[l] = 1-sum(af[1:(l-1)])
    else
        af[l]=1
    if (af[l]<=0)
    {
        af=rep(1,length(af))/length(af)
    } 
    rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@allele.freqs[[input$focalLoc]] <- af
})


observeEvent(input$changeNumAll,{
    rValues$ssClass@scenarios[[rValues$scenarioNumber]]@simulator.params@num.alleles <-
        rep(input$constNumAll,rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci)
})
