mig.mat <- reactive({
    if (is.null(input$distfun)) {dfun <- dexp} else {dfun <- get(input$distfun)}

    print(input$rows)
    print(input$cols)
    print(input$numpops)
    mat <- scenario.mig.matrix(
        h=input$numpops,
        h.dim=c(input$rows,input$cols),
        mig.model=input$migModel,
        distance.fun=dfun
        )
    mat$R.int * input$migRate
})


inmat <- reactive({
    mat <- matrix(0,input$numpops,input$numpops)
    for (row in 1:input$numpops)
        for (col in 1:input$numpops)
            {
                strng <- paste0("r",row,"c",col)
                print(input[[strng]])
                mat[row,col] <- input[[strng]]
            }
    mat
})

observe({
#    print(input$numpops)
#    print(paste("repopMig",input$repopulateMig))
    if (input$repopulateMig == 0) 
        return()
   isolate({
    output$migmat <-renderTable({

        mat <- mig.mat()
        retmat <- matrix("",dim(mat)[1],dim(mat)[2])
        for (row in 1:input$numpops)
            for (col in 1:input$numpops)
                {
                    intxt <- paste0("<input id='r",row,"c",col,"' class='input-tiny' type='number' value='",mat[row,col],"'>")
                    retmat[row,col] <-intxt
                }
        colnames(retmat) <- 1:input$numpops
        rownames(retmat) <- 1:input$numpops
        as.data.frame(retmat)
    }, sanitize.text.function = function(x) {x})
   })
})
