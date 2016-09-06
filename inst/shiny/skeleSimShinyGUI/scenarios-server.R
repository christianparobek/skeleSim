
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
##############################################
###############################################
mig.mat <- function(){
   if (debug()) print("in mig.mat")
    ##first check to see if there mig model is "user".  If so, and the matrix dimension has
    ##not changed, then multiply the matrix by migRate and return.  If the matrix dimension
    ##has changed, reset the matrix to an island model (but keep the scenario migration model
    ##equal to "user"

###set up some params (basically make short local names from long reactive object names)
    if (is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper))
        { #add in reasonable values for mig.helper
            rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper <-
                list(migModel="user",
                     migRate=1,
                     rows=as.integer(1),
                     cols=as.integer(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops),
                     distfun="dexp")
        } else {
            if (is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migModel))
                rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migModel <- "user"
            if (is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migRate))
                rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migRate <- 1
            if (is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$rows))
                rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$rows <- as.integer(1)
            if (is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$cols))
                rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$cols <-
                    as.integer(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops)
            if (is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$distfun))
                rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$distfun <- "dexp"
        }
    mmod <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migModel
    rws <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$rows
    cls <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$cols
    dfun <- get(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$distfun)
    numpop <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops
    
    if (is.null(mmod)) mmod <- "island"

    if (mmod=="user")
        {
            if (dim(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration[[rValues$migrationNumber+1]])[1]==numpop)
                { #use existing user matrix (proportional to migRate)
                    ret <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration[[rValues$migrationNumber+1]]
                } else { #make new user matrix
                    ret <- scenario.mig.matrix(
                        h=numpop,
                        mig.model="island")$R.int 
                }
        } else { #don't use a "user" migration model
            if  ((rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops)>1)
                {
                    if ((mmod %in% c("distance","twoD","twoDwDiagonal"))&(numpop!=rws*cls))#make sure 2d extent is correct
                        {
                            rValues$msg <- "Spatial model with non-rectangular extent: make sure rows*columns equals the number of pops"
                            ret <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration[[rValues$migrationNumber+1]]
                        } else if ((mmod %in% c("twoD","twoDwDiagonal"))&(min(rws,cls)<2))
                            {# not really a 2d landscape if all the pops are lined up....
                                rValues$msg <- "Spatial model 2D, but the landscape has either row or cols == 1"
                                ret <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration[[rValues$migrationNumber+1]]
                            } else { #looks good to proceed
                            rValues$msg <- NULL
                            mat <- scenario.mig.matrix(
                                h=numpop,
                                h.dim=c(rws,cls),
                                mig.model=mmod,
                                distance.fun=dfun
                                )
                            ret <- mat$R.int * rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mig.helper$migRate
                        }
                } else #when 1 pop, migration rate does not mean much
                    {
                        ret <- matrix(0,1,1)
                    }
        }
    if (debug()) print("about to return from mig.mat")
    ret
}


output$migmat <- renderUI({
    mat <- as.data.frame(mig.mat())
    matrixInput("migmat","Migration Matrix",
                isolate(mat))
})

pop.sizes <- reactive(
    {
        if (debug()) print("running popsize")
      if ((is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@pop.size))|
          (length(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@pop.size)!=
           rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops))
          ret <- rep(100,rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops)
      else 
          {
              ret <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@pop.size
          }
      matrix(ret,nrow=1)
    }
    )

output$popsize <- renderUI({
    if (debug()) print("creating popsize vector")
    psz <- as.data.frame(pop.sizes())
    matrixInput("psvec","Vector of population sizes",
                isolate(psz))
})

samp.sizes <- reactive(
    {
      if (debug()) print("running sampsize")
      if ((is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@sample.size))|
          (length(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@sample.size)!=
           rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops))
          ret <- rep(20,rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops)
      else 
          {
              ret <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@sample.size
          }
      matrix(ret,nrow=1)
    }
)

output$sampsize <- renderUI({
    if (debug()) print("creating popsize vector")
    ssz <- as.data.frame(samp.sizes())
    matrixInput("ssvec","Vector of sample sizes",
                isolate(ssz))
})

observeEvent(input$resetMutRate,{      
    if (length(grep("Gamma",input$specifyMutRate))>0)
    {
    rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mut.rate <-getGammaMutRates(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci,input$gammaMean,input$gammaStd)  
    } else {
     rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mut.rate <- rep(input$ConstMutRate,
                                                                         rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci)
    }
})

mut.rates <- reactive(
    {

      if (debug()) print("running sampsize")
      if ((is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mut.rate))|
          (length(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mut.rate)!=
           rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci))
          ret <- rep(1e-5,rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci)
      else 
          {
              ret <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mut.rate
          }
      matrix(ret,nrow=1)
    }
)

output$mutrate <- renderUI({
    mr <- as.data.frame(mut.rates())
    matrixInput("mutvec","Vector of mutation rates",
                isolate(mr))
})


output$scenario <- renderText({
    cat(str(rValues$ssClass@scenarios[[rValues$scenarioNumber]]))
})

