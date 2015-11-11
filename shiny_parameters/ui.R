#source("setup.R")

shinyUI(
  navbarPage(
    "skelesim",
      tabPanel("Intro questions",
               sidebarLayout(
                   sidebarPanel(
                                        #                 textInput("simname", "Simulation Name:", "Sim Parameters #1"),
                       checkboxInput("snps", label = "Do you have SNP data?", value = FALSE),
                       checkboxInput("non.diploid", label = "Is your data other than diploid?", value = FALSE),
                       checkboxInput("marker.num", label = "Do you want to simulate many markers?", value = FALSE),
                       checkboxInput("pop.size", label = "Do you have large population sizes?", value = FALSE),
                       checkboxInput("complex.hist", label = "Do you have a complex history to simulate?", value=FALSE),
                       checkboxInput("deep.time", label = "Are you looking at deep time frames", value = FALSE),
                       checkboxInput("demography", label = "Do you want to include demography?", value = FALSE),
                       checkboxInput("management", label = "Does your question involve management decisions?", value = FALSE),
                       checkboxInput("completion.time", label = "Do you need a short completion time", value = FALSE),
                       checkboxInput("computer", label = "Do you have large computer capacity?", value = FALSE),
                                        # for the file uploader
                       fileInput("file", label = h3("OR choose file to upload"))
                       ),
                   
                   mainPanel(
                       includeMarkdown("helpfiles/help-questions.md"),
                       h3(textOutput("simname", container = span)),
                       tableOutput("values")
                       )
                   )),
                                        #
      tabPanel("General Config.",
             sidebarLayout(
                 sidebarPanel(
                     actionButton("readss","Read SS object"),
                     textInput("title", "Title",
                               value = "Project title"),
                     checkboxInput("quiet", "Quiet?",
                                   value = FALSE),
                     checkboxInput("coalescent", "Coalescent simulator?",
                                   value = TRUE),
                     numericInput("reps", "Number of simulation reps",
                                  value = 1),
                     textInput("wd", "Working directory for simulation",
                               value = "wdTest")
                     ),
                 mainPanel()
                 )),  ##comma is critical between panels
                                        #     
      tabPanel("Scenario Config.",
               sidebarLayout(
                   sidebarPanel(
                       numericInput("scenarioNumber", "Which scenario",value=1),
                       textInput("numpopsTxt", "Number of Populations",
                                    value = "1"),
                       numericInput("numloci", "Number of loci",
                                    value = 1),
                       numericInput("mutrate", "Mutation Rate",
                                    value = 1e-4),
                       br(),
                       numericInput("migRate", "Migration rate multiplier (no effect for model 'user')",value=1),
                       selectInput("migModel", "Migration Model",
                                   choices=c("island","stepping.stone.linear",
                                       "stepping.stone.circular","twoD","twoDwDiagonal","distance","user"),
                                   selected=1),
                       conditionalPanel(condition = "input.migModel == 'distance' || input.migModel == 'twoD' || input.migModel == 'twoDwDiagonal'", numericInput("rows", "Rows in a grid-shaped landscape",1)),
                       conditionalPanel(condition = "input.migModel == 'distance' || input.migModel == 'twoD' || input.migModel == 'twoDwDiagonal'", numericInput("cols", "Cols in a grid-shaped landscape",1)),
                       conditionalPanel(condition = "input.migModel == 'distance'",
                                        selectInput("distfun", "Distance function (must be an R function)",c("dexp")))
                       )
                   , #don't forget comma
                   mainPanel(
                       tabsetPanel(
                           tabPanel("Migration matrix",
                                    textOutput("msg"), #potential error messages
                                    br(),
                                    includeMarkdown("helpfiles/help-migration.md"),
                                    uiOutput("migmat"),
                                    plotOutput("networkPlot"))
                           ))
                   ))
      
      ,  ##comma is critical between panels
      
      tabPanel("Specific simulator params",
               sidebarLayout(
                   sidebarPanel(
                      numericInput("specScenNumber","Scenario number",1),
                       conditionalPanel(
                           condition = "input.coalescent == true",
                           h4("Fastsimcoal parameters")),
                       conditionalPanel(
                            condition = "input.coalescent == true",
                            checkboxInput("infSiteModel", "Infinite site model",
                                          value = FALSE))),
                   mainPanel(
                       conditionalPanel(
                           condition = "input.coalescent == true",
                           tabsetPanel(
#                               tabPanel("Simcoal history specification",tableOutput("simhistTbl")),
                               tabPanel("Simcoal History",
                                        includeMarkdown("helpfiles/help-history.md"),
                                        plotOutput("simhistPlot",
                                                   click= "histplotClick",
                                                   dblclick = "histplotDblclick"),
                                        tableOutput("simhistTbl")
                               )
                           )
                       ))
               ))
      , #don't forget the comma

    tabPanel(
      "Run Simulator",
      uiOutput("btnRun"),
      textOutput("txtRunStatus")
    )
      ,

    tabPanel(
      "Current ssClass",
      tableOutput("ssClass")
    )
  )
)
