
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
### this one is a wrapper for the rmetasim function
### landscape.mig.matrix (remapped to scenario.mig.matrix)
### it has a lot of special-case handling built in.
###############################################
mig.mat <- function(){
   print("in mig.mat")
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
            if (dim(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration[[migrationNumber]])[1]==numpop)
                { #use existing user matrix (proportional to migRate)
                    ret <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration[[migrationNumber]]
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
                            ret <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration[[migrationNumber]]
                        } else if ((mmod %in% c("twoD","twoDwDiagonal"))&(min(rws,cls)<2))
                            {# not really a 2d landscape if all the pops are lined up....
                                rValues$msg <- "Spatial model 2D, but the landscape has either row or cols == 1"
                                ret <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration[[migrationNumber]]
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
    print("about to return from mig.mat")
    ret
}


output$migmat <- renderUI({
    matrixInput("migmat","Migration Matrix (please don't use +/- buttons at right [temporary])",
                as.data.frame(mig.mat()))
#                as.data.frame(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration))
})

pop.sizes <- function()
    {
        print("running popsize")
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

output$popsize <- renderUI({
    print("creating popsize vector")
    matrixInput("psvec","Vector of population sizes (please don't use +/- buttons at right [temporary])",
                as.data.frame(pop.sizes()))
})

samp.sizes <- function()
    {
        print("running sampsize")
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

output$sampsize <- renderUI({
    print("creating popsize vector")
    matrixInput("ssvec","Vector of sample sizes (please don't use +/- buttons at right [temporary])",
                as.data.frame(samp.sizes()))
})

mut.rates <- function()
    {
        print("running sampsize")
      if ((is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mut.rate))|
          (length(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mut.rate)!=
           rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci))
          ret <- rep(10e-5,rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.loci)
      else 
          {
              ret <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@mut.rate
          }
      matrix(ret,nrow=1)
    }

output$mutrate <- renderUI({
    matrixInput("mutvec","Vector of mutation rates (please don't use +/- buttons at right [temporary])",
                as.data.frame(mut.rates()))
#                as.data.frame(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration))
})

output$scenario <- renderText({
    cat(str(rValues$ssClass@scenarios[[rValues$scenarioNumber]]))
})

