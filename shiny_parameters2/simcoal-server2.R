mig.mat <- reactive({
    if (is.null(input$distfun)) {dfun <- dexp} else {dfun <- get(input$distfun)}
    mat <- scenario.mig.matrix(
        h=input$numpops,
        h.dim=c(input$rows,input$cols),
        mig.model=input$migModel,
        distance.fun=dfun
        )
    mat$R.int * input$migRate
})

inmat <- reactive({
    np <- input$numpops
    mat <- matrix(0,np,np)
    
    for (row in 1:np)
        for (col in 1:np)
            {
                strng <- paste0("r",row,"c",col)
                                        #            print(strng)
                                        #            print(input[[strng]])
                if (!is.null(input[[strng]]))
                    {
                        mat[row,col] <- input[[strng]]
                    } else {
                        mat[row,col] <- NA
                    }
            }
    mat
})


############################ observe expression, used for multiple purposes
observe({
    if (input$repopulateMig == 0) 
        return()
    isolate({
        output$migmat <-renderTable({
            mat <- mig.mat()
            retmat <- make.interactive.matrix(mat)
        }, sanitize.text.function = function(x) {x})
    })
})
