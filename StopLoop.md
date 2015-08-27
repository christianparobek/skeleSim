# Methods to stop simulation runs after a specified time
Michelle DePrenger-Levin  
Monday, April 06, 2015  

Want the ability to stop a simulation   
Will use RMetaSim to test


```r
library(rmetasim)
```

```
## Loading required package: ape
## Loading required package: ade4
## Loading required package: gtools
```

## Based on Time elapsed
Stop the repeat after the elapsed time exceeds a certain amount. Set to 1 but does not match 1 second...      
Using repeat{}


```r
# repeat and proc.time() ###
runs <- 1
ptm <- proc.time()
  repeat {
      exampleland <- landscape.new.example()
  is.landscape(exampleland)
  simex <- landscape.simulate(exampleland, 4)
  obshet <- landscape.obs.het(simex)
   # print(obshet)
  runs <- runs + 1
 newptm <- proc.time() - ptm
 # print(newptm)
if(newptm[3] >= 1) break()
  }
proc.time()
```

```
##    user  system elapsed 
##    1.91    0.21    3.10
```

```r
runs
```

```
## [1] 69
```

Learning a shiny to slide the time allowed would be nice. 

### This seems to measure actual time elapsed  
using repeat{} and break()


```r
start <- Sys.time()
repeat{
  exampleland <- landscape.new.example()
  is.landscape(exampleland)
  simex <- landscape.simulate(exampleland, 4)
  obshet <- landscape.obs.het(simex)
   # print(obshet)
  runs <- runs + 1
  end <- Sys.time()
  #Time elapsed
  timel <- as.numeric(end-start)
 
if(timel >= 2) break()
  }
proc.time()
```

```
##    user  system elapsed 
##    3.94    0.21    5.13
```

```r
runs
```

```
## [1] 220
```


### Try using while() loop


```r
runs <- 1
ptm <- proc.time()
while(newptm[3] < 5) {
  exampleland <- landscape.new.example()
  is.landscape(exampleland)
  simex <- landscape.simulate(exampleland, 4)
  obshet <- landscape.obs.het(simex)
   # print(obshet)
  runs <- runs + 1
 newptm <- proc.time() - ptm
}
runs
```

```
## [1] 390
```



### Stackoverflow 'tic/toc' like MATLAB   
But I can't get it to work. 



```r
tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self"))
{
   type <- match.arg(type)
   assign(".type", type, envir=baseenv())
   if(gcFirst) gc(FALSE)
   tic <- proc.time()[type]         
   assign(".tic", tic, envir=baseenv())
   invisible(tic)
}

toc <- function()
{
   type <- get(".type", envir=baseenv())
   toc <- proc.time()[type]
   tic <- get(".tic", envir=baseenv())
   print(toc - tic)
   invisible(toc)
}

runs <- 1
eltime <- 0
 tic()
while(toc() < 5){
exampleland <- landscape.new.example()
  is.landscape(exampleland)
  simex <- landscape.simulate(exampleland, 4)
  obshet <- landscape.obs.het(simex)
  runs <- runs + 1
}
```

```
## elapsed 
##       0
```

```r
runs
```

```
## [1] 1
```

