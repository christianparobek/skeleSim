This folder contains the shiny interface for skelesim.  As the codebase is structured now, everyting depends on this ui for 
organisation. In other words, there are no command-line approaches to skelesim at the moment.

To use,

1. From CRAN:
  1. make sure that the 'devtools' package is installed.
  2. make sure that the 'shiny' package is installed
  3. make sure that the 'rmetasim' package is installed
  4. make sure that the 'shinyFiles' package is installed
2. From github:
  1. load devtools package -- library(devtools)
  2. install_github("rstudio/shiny-incubator")
3. Set the R working directory to  the main skelesim directory:
  1. type library(shiny)
  2. type runApp("shiny_parameters")
