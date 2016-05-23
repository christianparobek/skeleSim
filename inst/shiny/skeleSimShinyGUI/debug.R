###
### server code to be included into 'server.R'
### this implements the printout of the ssClass object
###

output$ssClass <- renderTable({data.frame(item=(capture.output(str(rValues$ssClass))))})

