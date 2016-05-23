library(shiny)

# Define UI for dataset viewer application
shinyUI(fluidPage(
  
  # Application title
  titlePanel("skeleSim: The Population Simulation Builder"),
  
  # Sidebar with controls to provide a caption, select a dataset,
  # and specify the number of observations to view. Note that
  # changes made to the caption in the textInput control are
  # updated in the output area immediately as you type
  sidebarLayout(
    sidebarPanel(
      
      # for building your own sim
      textInput("simname", "Simulation Name:", "Population Simulation #1"),
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
      
      includeMarkdown("help.md"),
      h3(textOutput("simname", container = span)),
      tableOutput("values")
      
    )
  )
))