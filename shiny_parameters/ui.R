source("setup.R")

shinyUI(navbarPage("skelesim",

                   tabPanel("General Config.",
                            sidebarLayout(
                                sidebarPanel(
                                    textInput("title", "Title",
                                              value = "Project title"),

                                    textInput("date", "Date",
                                              value = Sys.time()),
                                    checkboxInput("quiet", "Quiet",
                                                  value = FALSE),
                                    checkboxInput("coalescent", "Coalescent Simulator?",
                                                  value = TRUE),
                                    numericInput("reps", "Number of simulation reps",
                                                 value = 1),
                                    numericInput("timing", "Timing",
                                                 value = 2),
                                    textInput("simfunc", "Simulation Function",
                                              value = "fsc.run"),
                                    textInput("wd", "Simulation working directory",
                                              value = "test.wd")
                                    ),
                                mainPanel()
                                )),
                   
                   tabPanel("Scenario Config.",
                            sidebarLayout(
                                sidebarPanel(
                                    numericInput("numpops", "Number of populations",
                                                 value = 3),
                                    numericInput("numloci", "Number of loci",
                                                 value = 1),
                                    numericInput("mut.rate", "Mutation Rate",
                                                 value = 1e-4),
                                    checkboxInput("repopulateMig","Rewrite migration matrix",TRUE),
                                    selectInput("migModel", "Migration Model",
                                                c("island","stepping.stone.linear",
                                                  "stepping.stone.circular","twoD","twoDwDiagonal","distance")),
                                    numericInput("migRate", "Migration rate",1),
                                    uiOutput("rows"),
                                    uiOutput("cols"),
                                    uiOutput("distanceFun")
                                    ),
                                mainPanel(
                                    tabsetPanel(
                                        tabPanel("Migration matrix",tableOutput("migmat")),
                                        tabPanel("Migration graph",plotOutput("networkPlot"))
                                        ))
                                )),
                   
                   tabPanel("Specific simulation config.",
                            sidebarLayout(
                                sidebarPanel(
                                    uiOutput("infsitesUI"),
                                    uiOutput("simhistUI")
                                    ),
                                
                                mainPanel(
                                    conditionalPanel(
                                        condition = "input.coalescent == true",
                                        tabsetPanel(
                                            tabPanel("Simcoal history specification",tableOutput("simhistTbl")),
                                            tabPanel("Graphical view",
                                                     plotOutput("simhistPlot",
                                                                click= "histplotClick",
                                                                dblclick = "histplotDblclick"))
                                            )
                                        )
                                    ))
                            )
                   )
        )
