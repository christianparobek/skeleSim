
#########################################
####### SOME TEST VISUALIZATIONS ########
#########################################

output$testViz1 <- renderPlot({
  input$stat1
  plot(1:10, 1:10, main = "First Statistic")
})

output$testViz2 <- renderPlot({
  input$newplot
  # Add a little noise to the cars data
  cars2 <- cars + rnorm(nrow(cars))
  plot(cars2)
})

#########################################
######### FOR GLOBAL STATISTICS #########
#########################################

## a list of some analyses we will want to graph
all_analyses <- c("Chi2", "Fst", "PHIst", "Tajimas.D", "heterozygosity", "Gst")

## given the current scenario, which of these analyses are actually present?
current_analyses <- all_analyses[all_analyses %in% colnames(result$params@analysis.results[[1]][[1]][,,1])]

## This is just for some practice for my own understanding
something <- current_analyses[1]
output$sometext <- renderPrint({current_analyses})


## setup the number of blocks that I'll plot in
output$plot_global <- renderUI({
  
  plot_output_list <- lapply(1:length(current_analyses), function(i) {
    plotname <- paste("plot", i, sep="")
    plotOutput(plotname, height = 280, width = 400)
  })
  
  # Convert the list to a tagList - this is necessary for the list of items
  # to display properly.
  do.call(tagList, plot_output_list)
})

for (i in 1:length(current_analyses)) {
  # Need local so that each item gets its own number. Without it, the value
  # of i in the renderPlot() will be the same across all instances, because
  # of when the expression is evaluated.
  local({
    my_i <- i
    plotname <- paste("plot", my_i, sep="")
    stat <- current_analyses[i]
    
    output[[plotname]] <- renderPlot({
      column <- which(colnames(result$params@analysis.results[[as.numeric(input$scenario)]][[1]][,,1]) == stat)
      data <- result$params@analysis.results[[as.numeric(input$scenario)]][[1]][1,column,]
      hist(data, main = stat)
    })
  })
}

