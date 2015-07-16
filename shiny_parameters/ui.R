source("setup.R")

shinyUI(
  navbarPage(
    "skelesim",

    tabPanel(
      "File",
      navlistPanel(
        tabPanel(
          "Load parameters",
          fileInput("fileParams", h4("Choose .Rdata File")),
          uiOutput("uiSelectParamObj"),
          textOutput("txtSelectedTitle")
        ),
        tabPanel(
          "Save parameters",
          #  !!!! GET 'label' FROM ssClass@label
          uiOutput("txtChangeTitle"),
          textInput("txtTitle", label = h4("Title"), value = ssClass@title),
          h4("Save parameter file"),
          uiOutput("uiBtnSaveParams"),
          textOutput("txtSaveStatus")
        )
      )
    ),

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

    tabPanel("General Config.",
             sidebarLayout(
               sidebarPanel(
                 uiOutput("titleUI"),
                 uiOutput("quietUI"),
                 uiOutput("coalescentUI"),
                 uiOutput("repsUI"),
                 uiOutput("wdUI")
               ),
               mainPanel()
             )),

    tabPanel("Scenario Config.",
             sidebarLayout(
               sidebarPanel(
                 uiOutput("scenarioNumberUI"),
                 uiOutput("numpopsUI"),
                 uiOutput("numlociUI"),
                 uiOutput("mutrateUI"),
                 br(),
                 actionButton("repopulateMig","Rewrite migration matrix"),
                 uiOutput("migmodelUI"),
                 uiOutput("migrateUI"),
                 uiOutput("rows"),
                 uiOutput("cols"),
                 uiOutput("distanceFun")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Migration matrix",
                            includeMarkdown("helpfiles/help-migration.md"),
                            tableOutput("migmat")),
                   tabPanel("Migration graph",
                            plotOutput("networkPlot"))
                   #                                        ,
                   #                                        tabPanel("debug",
                   #                                                 textOutput("scenDebug"))
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
                              includeMarkdown("helpfiles/help-history.md"),
                              plotOutput("simhistPlot",
                                         click= "histplotClick",
                                         dblclick = "histplotDblclick"))
                   )
                 )
               ))
    ),

    tabPanel(
      "Run Simulator",
      uiOutput("btnRun"),
      textOutput("txtRunStatus")
    )
  )
)
