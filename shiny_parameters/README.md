This folder contains the shiny interface for skeleSim.  As the codebase is structured now, everything depends on this ui for 
organisation. In other words, there are no command-line approaches to skelesim at the moment.

Right now (11/23/15) the app only supports  fastsimcoal simulations. Rmetasim simulations are in development.

To use,

1. From CRAN:
  1. make sure that the 'devtools' package is installed.
  2. make sure that the 'shiny' package is installed
  3. make sure that the 'rmetasim' package is installed
  4. make sure that the 'shinyFiles' package is installed
2. From github.com:
  1. load devtools package -- library(devtools)
  2. install_github("rstudio/shiny-incubator")
3. Set the R working directory to  the main skelesim directory:
  1. type library(shiny)
  2. type runApp("shiny_parameters")

**INSTALLING fastSimcoal2**

fastSimcoal2 binaries are availble at http://cmpg.unibe.ch/software/fastsimcoal2/

Make sure that you install the executable file in a place in your binary file search path.  This is critical.

**WINDOWS**

To use on windows, install the latest version of R and RStudio (obviously).  Then:

1. ensure that the fastsimcoal2 executable file is in your path (repeat of admonition above)
2. ensure that all of the R executables are in your path.
  -you can add to your path using the control panel in windows
  -alternatively you can download windows batchfiles that make this work more intuitively.  See https://code.google.com/p/batchfiles/ or http://cran.r-project.org/contrib/extra/batchfiles/ to find a copy of 'R.bat' if you put this file in your path, it handles searches for R on your windows disk.
3. then follow the instructions above.,,
