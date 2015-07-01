source("mig.matrix.R")
options(shiny.trace = F)  # cahnge to T for trace
require(shiny)
require(igraph)

shinyServer(function(input, output,session) {

                mig.mat <- reactive({
                    mat <- scenario.mig.matrix(
                        h=input$numpops,
                        mig.model=input$migModel
                        )
                    mat$R.int                        
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
                    if (input$repopulateMig == 0) 
                        return()
                    isolate({
                        output$migmat <-renderTable({
                            mat <- mig.mat()
                            retmat <- matrix("",dim(mat)[1],dim(mat)[2])
                            for (row in 1:input$numpops)
                                for (col in 1:input$numpops)
                                    {
#                                        intxt <- paste0("<input id='r",row,"c",col,"' class='shiny-bound-input' type='number' value='",mat[row,col],"'>")
                                        intxt <- paste0("<input id='r",row,"c",col,"' class='input-tiny' type='number' value='",mat[row,col],"'>")
#                                        print(intxt)
                                        retmat[row,col] <-intxt
                                }
                            colnames(retmat) <- 1:input$numpops
                            rownames(retmat) <- 1:input$numpops
                            as.data.frame(retmat)
                        }, sanitize.text.function = function(x) {x})
                    })
                })
                
                output$tbl <- renderTable({
                    inmat()
                    
                })
                
                output$txt <- renderText({
                                        #    names(input)
                                        #    c(input$coalescent,input[["coalescent"]])
                })

                
                output$distPlot <- renderPlot({
                    hist(rnorm(as.numeric(1000)))
                })

                output$networkPlot <- renderPlot({
                    mat <- as.matrix(inmat())
                    print(mat)
                    grph <-graph.adjacency(t(mat),weighted=T)
#                    autocurve.edges(grph)
                    plot(grph,
                         edge.label=round(E(grph)$weight, 2),
                         edge.curved=T)
                    
                })
})
