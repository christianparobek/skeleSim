output$testViz1 <- renderPlot({
  input$stat1
  plot(1:10, 1:10, main = "First Statistic")
})

output$testViz2 <- renderPlot({
  plot(1:10, 10:1, main = "Second Statistic")
})

output$plot <- renderPlot({
  input$newplot
  # Add a little noise to the cars data
  cars2 <- cars + rnorm(nrow(cars))
  plot(cars2)
})

output$fst <- renderPlot({
  data <- result$params@analysis.results[[as.numeric(input$scenario)]][[1]][1,3,]
  hist(data, main = "FST")
})

  
  
  
  