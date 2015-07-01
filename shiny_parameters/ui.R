library(shiny)
require(shinyIncubator)
                                        #library(shinysky)

shinyUI(fluidPage(

  titlePanel("Skelesim Configure"),

  sidebarLayout(
      sidebarPanel(
          radioButtons("modify","Edit the following parameters:",
                       c("Everything is good" = "dontModify",
                         "Simulation-level" = "paramModify",
                         "Scenario-level" = "scenarioModify",
                         "Simulator-specific" = "simulationModify")),

          br(), ##a little space 

          #panels for params = T
          conditionalPanel(
            condition = "input.modify == 'paramModify'",
              textInput("title", "Title",
                        value = "Simulation Title")
              ),
          conditionalPanel(
            condition = "input.modify == 'paramModify'",
              textInput("date", "Date",
                        value = Sys.time())
              ),
         conditionalPanel( 
            condition = "input.modify == 'paramModify'",
              checkboxInput("quiet", "Quiet",
                            value = FALSE)
              ),
         conditionalPanel( 
            condition = "input.modify == 'paramModify'",
              checkboxInput("coalescent", "Coalescent Simulator?",
                            value = TRUE)
              ),

          conditionalPanel( 
            condition = "input.modify == 'paramModify'",
              numericInput("reps", "Number of simulation reps",
                            value = 1)
              ),

          conditionalPanel( 
            condition = "input.modify == 'paramModify'",
              numericInput("timing", "Timing",
                            value = 2)
              ),

          conditionalPanel( 
            condition = "input.modify == 'paramModify'",
              textInput("simfunc", "Simulation Function",
                            value = "fsc.run")
              ),

          conditionalPanel( 
            condition = "input.modify == 'paramModify'",
              textInput("wd", "Simulation working directory",
                            value = "test.wd")
              ),


          

          #panels for scenario = T
          conditionalPanel(
            condition = "input.modify == 'scenarioModify'",
              numericInput("numpops", "Number of populations",
                        value = 3)
              ),
          conditionalPanel(
              condition = "input.modify == 'scenarioModify'",
              numericInput("numloci", "Number of loci",
                           value = 1)
              ),
          conditionalPanel(
              condition = "input.modify == 'scenarioModify'",
              numericInput("mut.rate", "Mutation Rate",
                           value = 1e-4)
              ),
          conditionalPanel(
              condition = "input.modify == 'scenarioModify'",
              selectInput("migModel", "Migration Model",
                           c("island","stepping.stone.linear",
                             "stepping.stone.circular"))
              ),
          conditionalPanel(
              condition = "input.modify == 'scenarioModify'",
              checkboxInput("repopulateMig", "start new matrix", TRUE)
              ),
          

          ##simulation specific parameters, dual tests sim edit and coalescent or not
          ## first set up the fastsimcoal parameters by checking if input.coalescent ==T
          conditionalPanel(
              condition = "input.modify=='simulationModify' & input.coalescent==true",
              checkboxInput("infSiteModel", "Infinite site model",
                            value = FALSE)
              )        
          ),

      mainPanel(
          textOutput("txt"),

          conditionalPanel(
              condition = "input.modify == 'scenarioModify'",
              h3("Migration matrix")
              ),

          conditionalPanel(
              condition = "input.modify == 'scenarioModify'",
              tableOutput("migmat")
              ),

          conditionalPanel(
              condition = "input.modify == 'scenarioModify'",
              h3("Graphical representation")
              ),
          conditionalPanel(
              condition = "input.modify == 'scenarioModify'",
              plotOutput("networkPlot")
              )
        )
      )
    ))

