
shinyUI(
  navbarPage(
    "skelesim",
    tabPanel(
        "Actions",
        h4("File operations" ),
        shinyFilesButton("fileParams","Read skelesim file","Select saved skelesim file", FALSE),
        uiOutput("uiSelectParamObj"),
        h5(textOutput("txtObjLabel")),
        shinySaveButton("ssClassSave","Save skelesim parameters to file","Save parameter file",filetype=list(ssClass=c("rdata","Rdata","rda"))),
        textOutput("txtSaveStatus"),
        h4("Run simulation"),
        shinyDirButton("workFolder","Select Simulation Root Directory","Set Simulation Root Directory",FALSE),
        br(),
        actionButton("btnRun","Run simulation"),
        br(),
        h4("Post simulation"),
        actionButton("quitbtn","Quit App")
    ),
    tabPanel("Intro questions",
             sidebarLayout(
                 sidebarPanel(
###                    textInput("simname", "Simulation Name:", "Sim Parameters #1"),
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
### Set parameters for the skeleSim class
    tabPanel("General Conf",
             sidebarLayout(
                 sidebarPanel(
                                        #                     actionButton("readss","Read SS object"),
                     textInput("title", "Title",
                               value = ""),
                     dateInput("date","Date"),
                     checkboxInput("quiet", "Quiet?",
                                   value = FALSE),
                     checkboxInput("coalescent", "Coalescent simulator?",
                                   value = FALSE),
                     numericInput("reps", "Number of simulation reps",
                                  value = 1),
                     checkboxGroupInput("analysesReq","Type of analyses for each rep",
                                        c("Global"="Global", "Pairwise"="Pairwise", "Locus"="Locus")),
                     textInput("wd", "Subdirectory for actual simulation",
                               value = "wdTest")
                     
                 ),
                 mainPanel(h3("Parameters automatically set:"),
                           h4(textOutput("simulator")),
                           h4(textOutput("simfunc")),
                           h4(textOutput("simpath"))
                           )
             )),  ##comma is critical between panels
                                        #     
    tabPanel("Scenario Conf",
             sidebarLayout(
                 sidebarPanel(
                     numericInput("scenarioNumber", "Which scenario",value=1),
                     textInput("numpopsTxt", "Number of Populations",
                               value = "1"),
                     numericInput("numloci", "Number of loci",
                                  value = 1),
                     selectInput("loctype","Type of locus",choices=c("sequence","microsatellite"),selected="sequence"),
                     conditionalPanel(condition = "input.loctype == 'sequence'",
                                      numericInput("seqlen","Sequence length",value=100)),
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
                         tabPanel("Population characteristics",
                                  uiOutput("popsize"),
                                  uiOutput("sampsize")
                                  ),
                         tabPanel("Among population migration",
                                  h4(textOutput("msg")), #potential error messages
                                  br(),
                                  includeMarkdown("helpfiles/help-migration.md"),
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
                    numericInput("specScenNumber","Scenario number",1),

##################### simcoal specific
                    conditionalPanel(
                        condition = "input.coalescent == true",
                        h4("Fastsimcoal parameters")),
                    conditionalPanel(
                        condition = "input.coalescent == true",
                        checkboxInput("infSiteModel", "Infinite site model",
                                      value = FALSE)),
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
                    )

                    
                ),
                mainPanel(
                    conditionalPanel(
                        condition = "input.coalescent == true",
                        tabsetPanel(
                            tabPanel("Simcoal History",
                                     includeMarkdown("helpfiles/help-history.md"),
                                     plotOutput("simhistPlot",
                                                click = "histplotClick",
                                                dblclick = "histplotDblclick"),
                                     tableOutput("simhistTbl")
                                     )
                           ,
                            tabPanel("Growth rates",
                                     uiOutput("growthrate")
                                     ),
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
                                     textOutput("leading"),
                                     tableOutput("lefkovitch"),
                                     matrixInUI("survmat"),
                                     matrixInUI("repmat"),
                                     matrixInUI("malemat")
                                     ),
                             tabPanel("rmetasim locus details",
                                     includeMarkdown("helpfiles/rmetasim-loci.md"),
                                     radioButtons("numfreq",NULL,c("use allele numbers","specify full allele frequency distribution")),
                                     vectorInUI("numall"),
                                     textOutput("afreqs")
                                     )

                        )
                    )
                    
                )
            )
            )
 , #don't forget the comma
   
   tabPanel(
       "Current ssClass",
       tableOutput("ssClass")
   )
  )
)
