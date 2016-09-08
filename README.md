[![Travis-CI Build Status](https://travis-ci.org/christianparobek/skeleSim.svg?branch=master)](https://travis-ci.org/christianparobek/skeleSim)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/christianparobek/skeleSim?branch=master&svg=true)](https://ci.appveyor.com/project/christianparobek/skeleSim)

# skeleSim

## Description

*skeleSim* is a tool to guide users in choosing appropriate simulations, setting parameters, calculating summary genetic statistics, and organizing data output, all within the R environment. *skeleSim* is designed to be an extensible environment that can 'wrap' around any simulation software to increase its accessibility and use.
    
## Installation

To install the latest version from GitHub:

```r
# make sure you have Rtools installed
if (!require('devtools')) install.packages('devtools')
# install from GitHub
devtools::install_github('christianparobek/skeleSim', build_vignettes = TRUE)
```

## Execution

To run the shiny app:

```r
# Load the skeleSim package
library(skeleSim)

# Run the app
skeleSimGUI()
```

## Contact

* submit suggestions and bug-reports: <https://github.com/christianparobek/skeleSim/issues>
* send a pull request: <https://github.com/christianparobek/skeleSim/>

## version 0.9.5

* First CRAN release.