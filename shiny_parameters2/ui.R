source("setup.R")

shinyUI(fluidPage(

  titlePanel("Skelesim Configure"),

  sidebarLayout(
      sidebarPanel(
          radioButtons("modify",label="Edit the following parameters:",
                       choices=c("Everything is good, do nothing" = "dontModify",
                           "Simulation-level" = "paramModify",
                           "Scenario-level" = "scenarioModify",
                           "Simulator-specific" = "simulationModify")),
          
          br(), ##a little space 
          
          uiOutput("titleui"),
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


          
          uiOutput("numpopPanel"),
          #panels for scenario = T
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
                             "stepping.stone.circular","twoD","twoDwDiagonal","distance"))
              ),
          uiOutput("rows"),
          uiOutput("cols"),
          uiOutput("migrationRate"),
          uiOutput("distanceFun"),
          conditionalPanel(
              condition = "input.modify == 'scenarioModify'",
              checkboxInput("repopulateMig", "start new matrix", TRUE)
              ),
      
          ##simulation specific parameters, dual tests sim edit and coalescent or not
          ## first set up the fastsimcoal parameters by checking if input.coalescent ==T

          uiOutput("infsitesUI"),
          uiOutput("simhistUI")

          ),

      mainPanel(
#          textOutput("txt"),
          conditionalPanel(
              condition = "input.modify == 'scenarioModify'",
              tabsetPanel(
                  tabPanel("Migration matrix",tableOutput("migmat")),
                  tabPanel("Migration graph",plotOutput("networkPlot"))
                  )
              ),
          conditionalPanel(
              condition = "input.modify == 'simulationModify'",
              tabsetPanel(
                  tabPanel("Simcoal history specification",tableOutput("simhistTbl")),
                  tabPanel("Graphical view",plotOutput("simhistPlot",
                                                       click= "histplotClick",
                                                       dblclick = "histplotDblclick"))
                  )
              )
        )
      )
    ))

