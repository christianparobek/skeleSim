##
##
##generic vector and  matrix io functions that can abstract how we do matrices in skelesim gui
##
##
vectorInUI.js <- function(id)
{
    ns <- NS(id)
    uiOutput(ns("vec"))    
}

vectorIn.js <- function(input,output,session,vec,label="Vector")
{
    output$vec <- renderUI({
        ns <- session$ns
        matrixInput(ns("vec"),label,
                    as.data.frame(matrix(vec,nrow=1)))
        
    })
  return(
      reactive({
          input$vec
      })
         )
}



vectorInRowUI <- function(id)
{
    ns <- NS(id)
    tableOutput(ns("vec"))    
}

vectorInRow <- function(input,output,session,vec,label="Vector")
{
    #produces a row style vector input
    
    output$vec <- renderTable({
        mat <- matrix(vec,nrow=1)        
        ns <- session$ns
        cols <- dim(mat)[2]
        retmat <- matrix("",nrow=1,ncol=cols)
        for (col in 1:cols)
                {
                    id <- ns(paste0("r",1,"c",col))
                    val <- as.character(mat[1,col])
                    if (is.na(val)) val <- ""
                    intxt <- paste0("<input id='",id,"' class='input-tiny' type='number' value='",
                                    val,"'>")
                    retmat[1,col] <-intxt
                }
        rownames(retmat) <- label
        colnames(retmat) <- 1:cols
        as.data.frame(retmat)
    }, sanitize.text.function = function(x) {x})
  return(
      reactive({
          ns <- session$ns
          if (!(is.null(req(vec))))
          {
              ret <- rep(0,length(vec))
              for (col in 1:length(vec))
              {
                  inputname <- paste0("r",1,"c",col)
                  # print(input[[ns(inputname)]])
                  ret[col] <- ifelse (is.null(input[[paste0("r",1,"c",col)]]),NA,input[[paste0("r",1,"c",col)]])
              }
              # print(ret)
              ret
          }
      })
  )
}
vectorInColUI <- function(id)
{
    ns <- NS(id)
    tableOutput(ns("vec"))    
}

vectorInCol <- function(input,output,session,vec,label="Vector")
{
    #produces a row style vector input
    
    output$vec <- renderTable({
        mat <- matrix(vec,nrow=1)
        # print(mat)
        ns <- session$ns
        rows <- dim(mat)[1]
        retmat <- matrix("",nrow=rows,ncol=1)
        for (row in 1:rows)
                {
                    
                    id <- ns(paste0("r",row,"c",1))
                    val <- as.character(mat[row,1])
                    if (is.na(val)) val <- ""
                    intxt <- paste0("<input id='",id,"' class='input-tiny' type='number' value='",
                                    val,"'>")
                    retmat[row,1] <-intxt
                }
        colnames(retmat) <- label
        rownames(retmat) <- 1:rows
        as.data.frame(retmat)
    }, sanitize.text.function = function(x) {x})
  return(
      reactive({
          ns <- session$ns
          if (!(is.null(req(vec))))
          {
              ret <- rep(0,length(vec))
              for (row in 1:length(vec))
              {
                  inputname <- paste0("r",row,"c",1)
                  # print(input[[ns(inputname)]])
                  ret[row] <- ifelse (is.null(input[[paste0("r",row,"c",1)]]),NA,input[[paste0("r",row,"c",1)]])
              }
              # print(ret)
              ret
          }
      })
  )
}


###########################################

vectorIn <- vectorIn.js
vectorInUI <- vectorInUI.js

##########################################
matrixInUI <- function(id)
{
    ns <- NS(id)
    uiOutput(ns("mat"))    
}

matrixIn <- function(input,output,session,mat,label="Matrix")
{
    #print("in matrixIn")
    output$mat <- renderUI({
        ns <- session$ns
        # print(mat)
        matrixInput(ns("mat"),label,as.data.frame(mat))
        
    })

    return(
        reactive({
#            print("in reactive from matrixIn")
            input$mat
        }))
    
}



#################
################# demography
################ 

    demomatUI <- function(id,label)
    {
        ns <- NS(id)
        tabPanel(label,
                 matrixInUI(ns("mat"))
                 )
    }

    demomat <- function(input,output,session,inmat,label="Survival")
    {
        # print("calling demomat")
        # print(paste("inmat",inmat))
        return(reactive({callModule(matrixIn,"mat",mat=inmat,label)()}))
    }


###################
