# AES 7/2/15
#
# these are renderUI statements that only make sense
# included into the skelesim parameter interface
# using source()

##########general parameter values
output$titleUI <- renderUI({
    val="Project title"
    if (!is.null(ssClass@title)) {val <- ssClass@title}
    textInput("title", "Title",
              value = "Project title")

})

output$quietUI <- renderUI({
    val=FALSE
    if (!is.null(ssClass@quiet)) {val <- ssClass@quiet}
    checkboxInput("quiet", "Quiet?",
              value = FALSE)

})

output$coalescentUI <- renderUI({
    val=TRUE
    if (!is.null(ssClass@simulator.type)) {val <- ifelse(ssClass@simulator.type=="c",T,F)}
    checkboxInput("coalescent", "Coalescent simulator?",
                  value = val)

})

output$repsUI <- renderUI({
    val=1
    if (!is.null(ssClass@num.reps)) {val <- ssClass@num.reps}
    numericInput("reps", "Number of simulation reps",
              value = val)

})


output$wdUI <- renderUI({
    val="wdTest"
    if (!is.null(ssClass@wd)) {val <- ssClass@wd}
    textInput("wd", "Working directory for simulation",
              value = val)

})




#scenario UIs  this allows the defaults to be set for particular scenarios

output$scenarioNumberUI <- renderUI({
    numericInput("scenarioNumber", "Which scenario",
             value = 1)
})

scenario.exists <- reactive({
    res <- FALSE
    if ((length(input$scenarioNumber)>0)&length(ssClass@scenarios)>0)
        if (input$scenarioNumber<=length(ssClass@scenarios))
            res <- TRUE
    res
})


output$numpopsUI <- renderUI({
    val <- 3
    if ((scenario.exists()) && (!is.null(ssClass@scenarios[[input$scenarioNumber]])))
        val <- ssClass@scenarios[[input$scenarioNumber]]@num.pops
    numericInput("numpops", "Populations",
             value = val)
})

output$numlociUI <- renderUI({
    val <- 3
    if ((scenario.exists()) && (!is.null(ssClass@scenarios[[input$scenarioNumber]])))
        {val <- ssClass@scenarios[[input$scenarioNumber]]@num.loci}
    numericInput("numloci", "Number of loci",
             value = val)
})



output$mutrateUI <- renderUI({
    val <- 1e-4
    if ((scenario.exists()) && (!is.null(ssClass@scenarios[[input$scenarioNumber]])))
        {
            val <- ssClass@scenarios[[input$scenarioNumber]]@mut.rate
        }
    numericInput("mut.rate", "Mutation Rate",
                 value = val)
})

output$migmodelUI <- renderUI({
    choices <- c("island","stepping.stone.linear",
                  "stepping.stone.circular","twoD","twoDwDiagonal","distance")
    val <- which(choices=="island")
    print(paste("input$scenarioNumber",input$scenarioNumber))
    print(paste("scenario.exists()",scenario.exists()))
#    browser()
    if ((scenario.exists()) && (!is.null(ssClass@scenarios[[input$scenarioNumber]])))
        {
          if (!is.null(ssClass@scenarios[[input$scenarioNumber]]))
              {
#                  browser()
                  val <- which(choices==ssClass@scenarios[[input$scenarioNumber]]@mig.helper$mig.model)
                  print(paste("mig model val",val))
              }
      }
    selectInput("migModel", "Migration Model",choices=choices,selected=val)
})

output$migrateUI <- renderUI(
    {
        val <- 1
    if ((scenario.exists()) && (!is.null(ssClass@scenarios[[input$scenarioNumber]])))
        if (!is.null(ssClass@scenarios[[input$scenarioNumber]]@mig.helper$mig.rate))
                {val <- ssClass@scenarios[[input$scenarioNumber]]@mig.helper$mig.rate}
        numericInput("migRate", "Migration rate",value=val)
    }
    )
                                        #
output$rows <- renderUI({
    if (!is.null(input$migModel))
        if (input$migModel%in%c("twoD","twoDwDiagonal","distance"))
            {
                numericInput("rows", "Rows in a grid-shaped landscape",2)
            }
})

output$cols <- renderUI({
    if (!is.null(input$migModel))
        if (input$migModel%in%c("twoD","twoDwDiagonal","distance"))
            {
                numericInput("cols", "Cols in a grid-shaped landscape",2)
            }
})


output$distanceFun <- renderUI({
    if (!is.null(input$migModel))
        if (input$migModel%in%c("distance"))
            {
                textInput("distfun", "Distance function (must be an R function)","dexp")
            }
})


################SIMCOAL SPECIFICS
output$infsitesUI <- renderUI({
    if (input$coalescent==T)
        {
            checkboxInput("infSiteModel", "Infinite site model",
                          value = FALSE)
        } else {
        }
})

output$simhistUI <- renderUI({
    if (input$coalescent==T)
        {
            checkboxInput("simhistflag", "Specify a demographic history for simulation?",
                          value = FALSE)
        } else {
        }
})


