                output$networkPlot <- renderPlot({
                    mat <- as.matrix(inmat())
                    Sys.sleep(0.05)                
                    grph <-graph.adjacency(t(mat),weighted=T)
                    if (input$migModel%in%c("twoD","twoDwDiagonal","distance"))
                        {
                            colfact=1.2
                            rowfact=1.2
                            layout=matrix(0,ncol=2,nrow=dim(mat)[1])
                            lycnt <- 1
                            for (col in 1:input$cols)
                                for (row in 1:input$rows)
                                    {
                                        layout[lycnt,] <- c(colfact * col,rowfact * row)
                                        lycnt <- lycnt+1
                                    }
                        } else {
                            layout=NULL
                        }
                            plot(grph,
                                 edge.label=round(E(grph)$weight, 2),
                                 edge.curved=T, layout=layout)
                    
                })
                
