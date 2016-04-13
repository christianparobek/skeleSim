
shinyUI(
  navbarPage(
    "skeleSim",
    tabPanel(
        "Actions",
        h4("File operations" ),
        shinyFilesButton("fileParams","Read skelesim file","Select saved skelesim file", FALSE),
        uiOutput("uiSelectParamObj"),
        h5(textOutput("txtObjLabel")),
        shinySaveButton("ssClassSave","Save skelesim parameters to file","Save parameter file",
                        filetype=list(ssClass=c("rdata","Rdata","rda"))),
        textOutput("txtSaveStatus"),
        h4("Run simulation"),
        shinyDirButton("workFolder","Set Simulation Root Directory (required before proceeding below)","Set Simulation Root Directory",FALSE),
        br(),br(),
    actionButton("btnSave","Save example inputs for each scenario"),
        br(),
    actionButton("btnRun","Run simulation"),
        br(),
        h4("After running a simulation..."),
        actionButton("quitbtn","Quit App")
    ),
    tabPanel("Intro questions",
             sidebarLayout(
                 sidebarPanel(
###                    textInput("simname", "Simulation Name:", "Sim Parameters #1"),
#                     checkboxInput("snps", label = "Do you have SNP data?", value = FALSE),
                     checkboxInput("non.diploid", label = "Is your data other than diploid?", value = FALSE),
                     checkboxInput("marker.num", label = "Do you want to simulate many markers?", value = FALSE),
                     checkboxInput("pop.size", label = "Do you have large population sizes?", value = FALSE),
                     checkboxInput("complex.hist", label = "Do you have a complex history to simulate?", value=FALSE),
                     checkboxInput("deep.time", label = "Are you looking at deep time frames", value = FALSE),
                     checkboxInput("demography", label = "Do you want to include demography?", value = FALSE),
                     checkboxInput("management", label = "Does your question involve management decisions?", value = FALSE),
                     checkboxInput("completion.time", label = "Do you need a short completion time", value = FALSE),
                     checkboxInput("computer", label = "Do you have large computer capacity?", value = FALSE)
#                    ,
                                        # for the file uploader
#                     fileInput("file", label = h3("OR choose file to upload"))
                 ),
                 
                 mainPanel(
                     includeMarkdown("helpfiles/help-questions.md"),
                     h3(textOutput("simname", container = span)),
                     tableOutput("values")
                 )
             )),
### Set parameters for the skeleSim class
    tabPanel("General Conf",
             sidebarLayout(
                 sidebarPanel(
                                        #                     actionButton("readss","Read SS object"),
                     textInput("title", "Title",
                               value = "title"),
                     dateInput("date","Date"),
                     checkboxInput("quiet", "Quiet?",
                                   value = FALSE),
                     checkboxInput("coalescent", "Coalescent simulator?",
                                   value = TRUE),
                     numericInput("reps", "Number of simulation reps",
                                  value = 1,min=1),
                     br(),
                     h5("Type of analyses requested for each replicate"),
                     checkboxGroupInput("analysesReq",label="",choices=c("Global"="Global", "Pairwise"="Pairwise", "Locus"="Locus")),
                     numericInput("NumPermReps","Number of permutations for significance tests during analysis",value=0,min=0),
                     br(),
                     textInput("wd", "Temporary subdirectory for simulation reps",
                               value = "wdTest"),
                     width=3 #width of general questions tab                    
                 ),
                 mainPanel(h3("SkeleSim general parameters:"),
                           includeMarkdown("helpfiles/skelesim-general-help.md")
#                           h4(textOutput("simulator")),
#                           h4(textOutput("simfunc")),
#                           h4(textOutput("simpath"))
                           )
             )),  ##comma is critical between panels
                                        #     
    tabPanel("Scenario Conf",
             sidebarLayout(
                 sidebarPanel(
                     numericInput("scenarioNumber", "Which scenario",value=1,min=1),
                     textInput("numpopsTxt", "Number of Populations",
                               value = "1"),
                     br(),
                     selectInput("migModel", "Migration Model",
                                 choices=c("island","stepping.stone.linear",
                                           "stepping.stone.circular","twoD","twoDwDiagonal","distance","user"),
                                 selected=1),
                     numericInput("migRate", "Migration rate multiplier (no effect for model 'user')",value=1),
                     conditionalPanel(condition = "input.migModel == 'distance' || input.migModel == 'twoD' || input.migModel == 'twoDwDiagonal'", numericInput("rows", "Rows in a grid-shaped landscape",value=1)),
                     conditionalPanel(condition = "input.migModel == 'distance' || input.migModel == 'twoD' || input.migModel == 'twoDwDiagonal'", numericInput("cols", "Cols in a grid-shaped landscape",value=1)),
                     conditionalPanel(condition = "input.migModel == 'distance'",
                                      selectInput("distfun", "Distance function (must be an R function)",choices=c("dexp"))),
                     selectInput("loctype","Type of locus",choices=c("sequence","microsatellite"),selected="sequence"),
                     conditionalPanel(condition = "input.loctype != 'sequence'", 
                                      numericInput("numloci", "Number of loci",
                                                   value = 1)),
                     
                     conditionalPanel(condition = "input.loctype == 'sequence'",
                                      numericInput("seqlen","Sequence length",value=100)),

                     width=3 #number between 1-12 for sidebar width on scenarios tab                     
                 )
               , #don't forget comma
                 mainPanel(
                     tabsetPanel(
                         tabPanel("Population characteristics",
                                  uiOutput("popsize"),
                                  uiOutput("sampsize")
                                  ),
                         tabPanel("Among population migration",
                                  
                                  includeMarkdown("helpfiles/help-migration.md"),
                                  br(),
                                  h4(textOutput("msg")), #potential error messages
                                  fluidRow(column(3,numericInput("mignum","Migration matrix number",min=0,value=0,width='100px')),
                                           column(4,textOutput("numMigMats"))),
                                  uiOutput("migmat"),
                                  plotOutput("networkPlot")
                                  ),
                         tabPanel("Locus characteristics",
                                  uiOutput("mutrate")
                                  )
                     ))
             ))
    
 ,  ##comma is critical between panels

   tabPanel(textOutput("simtext"),
            sidebarLayout(
                sidebarPanel(
                    numericInput("specScenNumber","Scenario number",value=1,min=1),

##################### simcoal specific
                    conditionalPanel(
                        condition = "input.coalescent == true",
                        h4("Fastsimcoal parameters")),
#                    conditionalPanel(
#                        condition = "input.coalescent == true",
#                        checkboxInput("infSiteModel", "Infinite site model",
#                                      value = FALSE)),
                    conditionalPanel(
                        condition = "input.coalescent == true",
                        uiOutput("simexec")
                    )
#################### end simcoal specific
                    ,
##################### rmetasim specific
                     conditionalPanel(
                        condition = "input.coalescent != true",
                        numericInput("self", "selfing rate",
                                     value = 0)  ,
                        checkboxInput("changedemo","Change the number of stages within populations (drastic)",
                                      value=FALSE)
                    ),
                    conditionalPanel(
                        condition = "input.coalescent != true",
                        conditionalPanel(condition = "input.changedemo == true",
                            numericInput("stages", "demographic stages per location",
                                     value = 2)
                            )
                    ),
                    conditionalPanel(
                        condition = "input.coalescent != true",
                        numericInput("gens", "generations to simulate",
                                     value = 50)
                    ),
                    width=3 #width 1-12 of simulator specific tab                    
                ),
                mainPanel(
                    conditionalPanel(
                        condition = "input.coalescent == true",
                        tabsetPanel(
                            tabPanel("Simcoal History",
                           #          includeMarkdown("helpfiles/help-history.md"),
#                                     plotOutput("simhistPlot",
#                                                click = "histplotClick",
#                                                dblclick = "histplotDblClick"),
                                     plotOutput("simhistPlot"),
                                     uiOutput("simhistEditTbl")
                                    ,
                                     actionButton("addHistEvent","Add a new historical event"),
                                     actionButton("removeLastEvent","Remove the last row")
                                     )
                           ,
                            tabPanel("Growth rates",
                                     uiOutput("growthrate")
                                     )
                           ,
                            tabPanel("Sample times",
                                     uiOutput("samptime")
                                     )
                        )
                    ),
                    conditionalPanel(
                        condition = "input.coalescent != true",
                        tabsetPanel(
                            tabPanel("Within-population Lefkovitch Matrix",
                                     includeMarkdown("helpfiles/lefkovitch.md"),
                                     tableOutput("lefkovitch"),
                                     h5(textOutput("leading")),
                                     br(),
                                     h4("Edit elements below to change matrix above"),
                                     matrixInUI("survmat"),
                                     matrixInUI("repmat"),
                                     matrixInUI("malemat")
                                     ),
                             tabPanel("rmetasim locus details",
                                     includeMarkdown("helpfiles/rmetasim-loci.md"),
                                     radioButtons("numfreq",NULL,c("numbers of alleles")),
                                     vectorInUI("numall")
#                                    ,
#                                     textOutput("afreqs")
                                     )

                        )
                    )
                    
                )
            )
            ) ,
   tabPanel(
       "Current ssClass",
       tableOutput("ssClass")
   )
# , #don't forget the comma


      ##### MAKE THIS TAB REACTIVE SOMEHOW... EITHER NOT SHOW BEFORE THERE ARE RESULTS, OR SAY "NOTHING TO SEE HERE"
#     tabPanel("Visualize",
#          sidebarLayout(
#               sidebarPanel(
#                 selectInput("scenario", label = h3("Choose Scenario to Visualize:"),
#                             choice = c("Scenario #1" = 1,
#                                        "Scenario #2" = 2,
#                                        "Scenario #3" = 3)),
#                             ### NEED TO MAKE THIS REACTIVE TO NUMBER OF SCENARIOS THERE ARE; CAN HAVE > 3
#                 actionButton("newplot", "New plot")),
#               mainPanel(
#                 tabsetPanel(
#                   tabPanel("Global Statistics",
#                            uiOutput("plot_global"),
#                            fluidRow(column(3, verbatimTextOutput("sometext")))
#                            ),
#                   tabPanel("Locus Statistics",
#                            plotOutput("testViz1")),
#                   tabPanel("Pairwise Statistics",
#                            plotOutput("testViz2"))))
#       )
    )
  )

