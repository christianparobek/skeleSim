mig.mat <- reactive({
    if (is.null(input$distfun)) {dfun <- dexp} else {dfun <- get(input$distfun)}
   if  (!is.null(input$numpops))
       {
           mat <- scenario.mig.matrix(
               h=input$numpops,
               h.dim=c(input$rows,input$cols),
               mig.model=ifelse(is.null(input$migModel),"island",input$migModel),
               distance.fun=dfun
               )
           mat$R.int * input$migRate
       }
})


inmat <- reactive({
    mat <- matrix(0,input$numpops,input$numpops)
    for (row in 1:input$numpops)
        for (col in 1:input$numpops)
            {
                strng <- paste0("r",row,"c",col)
                mat[row,col] <- ifelse(is.null(input[[strng]]),NA,input[[strng]])
            }
    mat
})

#observe({
#    print(ssClass@scenarios[[input$scenarioNumber]]@migration[[1]])
#    print(mig.mat())
#    ssClass@scenarios[[input$scenarioNumber]]@migration[[1]] <- mig.mat()
#    print(ssClass@scenarios[[input$scenarioNumber]]@migration[[1]])
#},priority=10)

####the intent here is for every change to a mig mat to be written to the global ssClass object
update.ssClass.mig <- observe({  
    if (scenario.exists())
        {
#            mat <- mig.mat()
            print("theres a scenario")
            print(paste("scenario number ",input$scenarioNumber))
            print(ssClass@scenarios[[input$scenarioNumber]]@migration[[1]])
            print(input$migModel)
            print(ssClass@scenarios[[input$scenarioNumber]]@mig.helper$mig.model)
            if ((sum(is.na(ssClass@scenarios[[input$scenarioNumber]]@migration[[1]]))>0)
                |(input$migModel!=ssClass@scenarios[[input$scenarioNumber]]@mig.helper$mig.model))
                {
                    print("writing matrix")
                    ssClass@scenarios[[input$scenarioNumber]]@migration[[1]] <<- mig.mat()
                } else {
                    print ("matrix not altered")
                }
            print("this is the current state of the matrix")
            print(ssClass@scenarios[[input$scenarioNumber]]@migration[[1]])
        }
})


htmlmat <- reactive({
        input$migModel
        if (!is.null(input$numpops))
            {
                mat <- ssClass@scenarios[[input$scenarioNumber]]@migration[[1]]
                print("inside of htmlmat")
                print(mat)
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
            }
    })

output$migmat <-renderTable({
    htmlmat()
}, sanitize.text.function = function(x) {x})

