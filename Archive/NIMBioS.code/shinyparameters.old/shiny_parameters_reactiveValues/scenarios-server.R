
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
mig.mat <- reactive({
print("in mig.mat")
    if (is.null(input$distfun)) {dfun <- dexp} else {dfun <- get(input$distfun)}
    if  (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops))
        {
            mat <- scenario.mig.matrix(
                h=rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops,
                h.dim=c(input$rows,input$cols),
                mig.model=ifelse(is.null(input$migModel),"island",input$migModel),
                distance.fun=dfun
                )
            ret <- mat$R.int * input$migRate
        }
print("about to return from mig.mat")
ret
})


inmat <- reactive({

    print("entering inmat")
print(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops)
    if(is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops))
        {
            print("null num.pops")
            mat <- matrix(1,1,1)
        } else {
            print("num.pops non-null")
            mat <- matrix(0,rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops,
                          rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops)
        }
    for (row in 1:rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops)
        for (col in 1:rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops)
            {
                strng <- paste0("r",row,"c",col)
                print(strng)
                val <- ifelse(is.null(input[[strng]]),NA,input[[strng]])
                mat[row,col] <- val
            }
    mat
})




####the intent here is for every change to a mig mat to be written to the global ssClass object
ssClass.scenario.mig <- reactive({  
    if (scenario.exists())
        {
            if (rValues$ssClass@scenarios[[input$scenarioNumber]]@mig.helper$migModel!="user")
                {
#                    browser()
                    mat <- mig.mat()
                }
            else if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration[[1]]))
                {
                    mat <- rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration[[1]]
#                    if (dim(mat)[1]!=input$numpops)
#                        mat <- matrix(0,input$numpops,input$numpops)
                } else {
                    mat <- inmat()
                }
            if (sum(is.na(mat))>0)
                {
                    mat <- matrix(0,rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops,rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops)
                } 
        } else
            {
                print("scenario not exist, making mig mat")
                mat <- mig.mat()
            }
    mat
})


htmlmat <- reactive({
        if (!is.null(rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops))
            {
                mat <- ssClass.scenario.mig()
                ##this line updates the rValues$ssClass !!!!
                rValues$ssClass@scenarios[[rValues$scenarioNumber]]@migration[[1]] <- mat

                retmat <- matrix("",dim(mat)[1],dim(mat)[2])
                for (row in 1:rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops)
                    for (col in 1:rValues$ssClass@scenarios[[rValues$scenarioNumber]]@num.pops)
                        {
                            intxt <- paste0("<input id='r",row,"c",col,"' class='input-tiny' type='number' value='",mat[row,col],"'>")
                            retmat[row,col] <-intxt
                        }
                colnames(retmat) <- 1:dim(retmat)[2]
                rownames(retmat) <- 1:dim(retmat)[1]
                as.data.frame(retmat)
            }
    })

output$migmat <-renderTable({
#    imat <- inmat()
    hmat <- htmlmat()
    hmat
}, sanitize.text.function = function(x) {x})

output$scenario <- renderText({
    cat(str(rValues$ssClass@scenarios[[rValues$scenarioNumber]]))
})

