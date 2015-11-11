
###############################################
#### User interface renderers
#scenario UIs  this allows the defaults to be set for particular scenarios


##########################################################################
scenario.exists <- reactive({
    res <- FALSE
    if (length(rValues$ssClass@scenarios)>0)
        if (rValues$scenarioNumber<=length(rValues$ssClass@scenarios))
            res <- TRUE
    res
})


###############################################
### functions that operate only in the server
mig.mat <- function(){
    print("in mig.mat")
    ##first check to see if there mig model is "user".  If so, and the matrix dimension has
    ##not changed, then multiply the matrix by migRate and return.  If the matrix dimension
    ##has changed, reset the matrix to an island model (but keep the scenario migration model
    ##equal to "user"

###set up some params (basically make short local names from long reactive object names)
    mmod <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migModel
    rws <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$rows
    cls <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$cols
    dfun <- get(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$distfun)
    numpop <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops
    
    if (is.null(mmod)) mmod <- "island"

    if (mmod=="user")
        {
            if (dim(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration[[1]])[1]==numpop)
                { #use existing user matrix (proportional to migRate)
                    ret <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration[[1]]
                } else { #make new user matrix
                    ret <- scenario.mig.matrix(
                        h=numpop,
                        mig.model="island")$R.int 
                }
        } else { #don't use a "user" migration model
            if  ((rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops)>1)
                {
                    mat <- scenario.mig.matrix(
                        h=numpop,
                        h.dim=c(rws,cls),
                        mig.model=mmod,
                        distance.fun=dfun
                        )
                    ret <- mat$R.int * rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migRate
                } else #when 1 pop, migration rate does not mean much
                    {
                        ret <- matrix(0,1,1)
                    }
        }
    print ("new ret")
    print(ret)
    print("about to return from mig.mat")
    ret
}


output$migmat <- renderUI({
    matrixInput("migmat","Migration Matrix (please don't use buttons at right [temporary])",
                as.data.frame(mig.mat()))
#                as.data.frame(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration))
})

output$scenario <- renderText({
    cat(str(rValues$ssClass@scenarios[[rValues$scenarioNumber]]))
})

