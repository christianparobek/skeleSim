
shinyUI(
  navbarPage(
    "skeleSim",
    tabPanel(
        "Actions",
        checkboxInput("ActionHelp","Check to show verbose help for this page",TRUE),
        h4("File operations" ),
        conditionalPanel(condition="input.ActionHelp",
                         includeMarkdown("helpfiles/file-ops.md")),
        shinyFilesButton("fileParams","Read skelesim object from file","Select saved skelesim file", FALSE),
        uiOutput("uiSelectParamObj"),
        h5(textOutput("txtObjLabel")),
        shinySaveButton("ssClassSave","Save skelesim object to file","Save parameter file",
                        filetype=list(ssClass=c("rdata","Rdata","rda"))),
        textOutput("txtSaveStatus"),
        h4("Run simulation"),
        conditionalPanel(condition="input.ActionHelp",
                         includeMarkdown("helpfiles/simulations.md")),
        shinyDirButton("workFolder","Set Simulation Root Directory (required before proceeding below)","Set Simulation Root Directory",FALSE),
        textOutput("rootText"),
        br(),
    actionButton("btnSave","Save example inputs for each scenario"),
        br(),
    actionButton("btnRun","Run simulation"),
    br(),
    h4("After submitting a simulation,"),
    h4("re-edit and submit/save another, or..."),
#    actionButton("resetButton","Reset skeleSim object to default"),
    actionButton("quitbtn","Quit skeleSimGUI()")
    ),
    tabPanel("Help Choosing Simulator",
             sidebarLayout(
                 sidebarPanel(
###                    textInput("simname", "Simulation Name:", "Sim Parameters #1"),
#                     checkboxInput("snps", label = "Do you have SNP data?", value = FALSE),
                     checkboxInput("non.diploid", label = "Is your data other than diploid?", value = FALSE),
                     checkboxInput("marker.num", label = "Do you want to simulate many markers?", value = FALSE),
                     checkboxInput("pop.size", label = "Do you have large population sizes?", value = FALSE),
                     checkboxInput("complex.hist", label = "Do you have a complex history to simulate?", value=FALSE),
                     checkboxInput("deep.time", label = "Are you looking at deep time frames", value = FALSE),
                     checkboxInput("demography", label = "Include within population demography?", value = FALSE),
                     checkboxInput("management", label = "Does your question involve management decisions?", value = FALSE),
                     checkboxInput("completion.time", label = "Do you need a short completion time", value = FALSE),
                     checkboxInput("computer", label = "Do you have large computer capacity?", value = FALSE)
#                    ,
                                        # for the file uploader
#                     fileInput("file", label = h3("OR choose file to upload"))
                 ),
                 
                 mainPanel(
                     includeMarkdown("helpfiles/help-questions.md"),
                     tableOutput("values"),
                     br(),
                     h4("Relative weighting to coalescent versus forward-time simulator"),
                     textOutput("coalscore"),
                     textOutput("forwardscore")
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
                     selectInput("scenarioNumber", "Which scenario",choices=1,selected=1),
                     actionButton("addScenario","Add a new scenario"),
                     br(),br(),
                     textInput("numpopsTxt", "Number of Populations",
                               value = "2"),
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
                     selectInput("loctype","Type of locus",choices=c("sequence","microsatellite","SNP"),selected="microsatellite"),

                     conditionalPanel(condition = "input.loctype == 'sequence'",
                                      numericInput("seqlen","Sequence length",value=100)),

                     textInput("numloci", "Number of loci",
                               value = "1"),

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
                                  selectInput("specifyMutRate","Populate mutation rate vector with:",selected=1,choices=c("Constant rate","Rates chosen at random from Gamma Distribution")),
                                  conditionalPanel(condition = "input.specifyMutRate == 'Constant rate'",
                                                   
                                  numericInput("ConstMutRate",
                                               "Constant mutation rate across loci",1e-04)),
                                  conditionalPanel(condition = "input.specifyMutRate != 'Constant rate'",
                                                   
                                  numericInput("gammaMean",
                                               "Mean of Gamma Distribution",1e-04)),
                                  conditionalPanel(condition = "input.specifyMutRate != 'Constant rate'",
                                                   
                                  numericInput("gammaStd",
                                               "StdDev of Gamma Distribution",1e-04)),
                                  actionButton("resetMutRate","RePopulate Mutation Rate Vector"),
                                  uiOutput("mutrate")
                                  )
                     ))
             ))
    
 ,  ##comma is critical between panels

   tabPanel(textOutput("simtext"),
            sidebarLayout(
                sidebarPanel(
                    selectInput("specScenNumber","Scenario number",choices=c(1),selected="1"),

##################### simcoal specific
                    conditionalPanel(
                        condition = "input.coalescent == true",
                        h4("Fastsimcoal parameters")),
#                    conditionalPanel(
#                        condition = "input.coalescent == true",
#                        {
#                            checkboxInput("coalesceStar", "Make all pops coalesce at max time in past",
#                                          value = FALSE)
#                        }
                    conditionalPanel(
                        condition = "input.coalescent == true",
#                        uiOutput("simexec")
                        selectInput("fscexec","Select fastsimcoal executable",
                                          choices="No compatible executable found in path")
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
                        numericInput("gens", "time-periods (years?) to simulate",
                                     value = 50,min=1)
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
#                                     radioButtons("numfreq",NULL,c("numbers of alleles")),
                                     numericInput("constNumAll","How many alleles per locus?",1,min=1),
                                     actionButton("changeNumAll","Replace number of alleles across all loci to number above"),
                                     helpText('The following vector gives the number of alleles at each of the loci specified',
                                              'in the "Scenario Conf" tab'),
                                     uiOutput("numall"),
                                     includeMarkdown("helpfiles/rmetasim-alleles.md"),
                                     uiOutput("focalLoc"),
                                     uiOutput("afreqLoc")

#                                    ,
#                                     textOutput("afreqs")
                                     )

                        )
                    )
                    
                )
            )
            ) 
 , #don't forget the comma

   tabPanel("Results",
            sidebarLayout(
                sidebarPanel(
                    radioButtons("scompare","Type of plots",c("Summarize Scenario"="s","Compare Scenarios"="c"),
                                       selected="s"),
                    conditionalPanel(
                        condition = "input.scompare == 's'",
                        numericInput("vizScenario","Scenario to plot",value=1,min=1,max=1)
                    ),

                    selectizeInput("gstatsel","Global statistics",selected=NULL,multiple=TRUE,choices="None (load completed simulation)"),
                    selectizeInput("lstatsel","Locus-level statistics",selected=NULL,multiple=TRUE,choices="None (load completed simulation)"),
                    selectizeInput("pstatsel","Pairwise statistics",selected=NULL,multiple=TRUE,choices="None (load completed simulation)")
                )
               ,
                mainPanel(
                    tabsetPanel(
                        tabPanel("Global Statistics",
                                 downloadButton("dlGlobal","Download .csv of global analyses"),
                                 plotOutput("globalMainPlot",height=500,width=500)
                                 ),
                        tabPanel("Locus Statistics",
                                 downloadButton("dlLocus","Download .csv of locus analyses"),
                                 plotOutput("locusMainPlot",height=900,width=600)
                                 ),
                        tabPanel("Pairwise Statistics",
                                 downloadButton("dlPairwise","Download .csv of pairwise analyses"),
                                 plotOutput("pairwiseMainPlot",height=600,width=900)
                                 )
                   )))),
   tabPanel(
       "Current ssClass",
       tableOutput("ssClass")
   )
  )
)

